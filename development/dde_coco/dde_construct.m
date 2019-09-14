function prob=dde_construct(prob,type,varargin)
%% constructor for dde-biftool problems
%% Inputs
%
% * prob: problem structure created by coco_prob()
% * type: a DDE-Biftool point type/kind. Permitted/tested are 'stst',
% 'hopf', 'fold', 'psol', 'POfold', 'torus', 'PD'.
% * 'id' (default ''): coco id used as prefix for all parameters, functions
% * 'funcs': functions structure containing r.h.s. (as created by
% set_funcs...)
% * 'info': a single stst point (for initial 'stst'), or a DDE-biftool
% branch structure with method, parameter, point (and optionally) tangent
% fields.Then the first element in the point field will be used as initial
% point
% * pnames: cell array or index struct for parameter names. If the
% parameter field in point is longer, the remaining names will be parnum
% where num are integers (e.g. par5, par6). option parbasename can
% overwrite this.
% * tangent ([]): if logical & true, use tangent information from info, if
% point structure, use this as tangent for coco_add_func
% * 'data': data structure from previous run, which replaces all defaults
%% Output
% updated problem structure
%
%$Id: dde_construct.m 346 2019-05-13 05:41:50Z jansieber $
%% optional inputs
default={'id','','funcs',[],'info',[],'pnames',{},'parbasename','par',...
    'tangent',[],'bd_point',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on','data');
%% construction of initial solution
% possibly redefine the funcs structure, create initial solution and
% parameters stored in typical "branch" structure of DDE-Biftool branch is
% supposed to contain one point and (optional) one tangent. This function
% calls, if present, dde_kind2type, where type is input, kind is field in
% options.info.point (or options.info). If not present it calls SetupType
% (Type capitalized) from ddebiftool_utilities.
data=setup_br(type,options.funcs,options.info,pass_on{:});
%% fix: 
% some nonlinear problems use non-square Jacobians in DDE-Biftool. We avoid
% this for coco.
if strcmp(data.info.method.point.preprocess,'dde_jac2square_preprocess')
    data.info.method.point.preprocess='';
end
if isempty(options.id)
    fid='dde';
else
    fid=options.id;
end
id=options.id;
data.id=id;
data.userfuncs=options.funcs;
data.ipar=1:length(data.info.point.parameter);
ind=dde_ind_from_point(data.info.point,data.ipar);
u0=dde_x_from_point(data.info.point,data.ipar);
t_info={};
if ~isempty(options.tangent) && islogical(options.tangent) && options.tangent
    utan=dde_x_from_point(data.info.tangent,data.ipar);
    t_info={'t0',utan};
elseif ~isempty(options.tangent) && isstruct(options.tangent)
    utan=dde_x_from_point(options.tangent,data.ipar);
    t_info={'t0',utan};    
elseif isfield(data.info,'tangent')
    data.info=rmfield(data.info,'tangent');
end
data.pnames=fill_pnames(options.pnames,options.parbasename,data.ipar);
pnames_global=coco_get_id(id,data.pnames);
data.userpars=setdiff(data.ipar,data.info.parameter.free);
data.extpars=data.info.parameter.free;
data=coco_func_data(data);
prob=coco_add_func(prob,fid,@dde_rhs,data,'zero','u0',u0,t_info{:}, 'remesh', @dde_remesh,'f+df');
uidx=coco_get_func_data(prob,fid,'uidx');
xnorm_name=coco_get_id(id,'xnorm');
prob=coco_add_func(prob,xnorm_name,@xnorm,struct('id',fid),'regular',xnorm_name,'uidx',uidx);
%% add all parameters of point as inactive to problem
if ~isempty(data.userpars)
    prob=coco_add_pars(prob,'pars',uidx(ind.parameter(data.userpars)),...
        pnames_global(data.userpars),'inactive');
end
if ~isempty(data.extpars)
    prob=coco_add_pars(prob,'extpars',uidx(ind.parameter(data.extpars)),...
        pnames_global(data.extpars),'active');
end
%% Reference point, used by some nonlinear problems (phase conditions etc)
prob = coco_add_slot(prob, fid, @dde_point_update, data, 'update'); 
%% save results, if requested put complete solutions into column of bd
prob = coco_add_slot(prob, fid, @coco_save_data, data, 'save_full');
if options.bd_point
    ptoutput=coco_get_id(id,'out');
    prob = coco_add_slot(prob, ptoutput, @bddat, data, 'bddat');
end
end
%%
function data=setup_br(type,funcs,info,varargin)
utype=[upper(type(1)),type(2:end)];
if strcmp(utype,'Stst')
    info_ind={'point',info};
else
    info_ind={info,1};
    constructor=['dde_',info.point.kind,'2',type];
    if exist(constructor,'file')
        data=dde_apply({'dde_',info.point.kind,['2',type]},funcs,info,varargin{:});
        return
    end
end
[out1,out2]=feval(['Setup',utype],funcs,info_ind{:},...
    'correc',false,'correction',false,varargin{:});
if ~isstruct(out2)
    data.info=out1;
    data.funcs=funcs;
else
    data.info=out2;
    data.funcs=out1;
end
end
%% monitoring function for solution norm
function [data, y] = xnorm(prob, data, u)
ptdata=coco_get_func_data(prob,data.id,'data');
pt=dde_point_from_x(u,ptdata.info.point(1),ptdata.ipar);
if isfield(pt,'x')
    y=norm(pt.x,2);
elseif isfield(pt,'profile')
    y=sqrt(p_dot(pt,pt));
end
end
%%
function pnames=fill_pnames(names,base,ind)
pnames=arrayfun(@(i)[base,num2str(i)],ind,'uniformoutput',false);
if isstruct(names)
    fnames=fieldnames(names);
    cind=struct2cell(names);
    [~,ism]=ismember([cind{:}],ind);
    pnames(ism)=fnames;
else
    pnames(1:length(names))=names;
end
%pnames=coco_get_id(parid,pnames);
end
%%
function [prob, stat, xtr] = dde_remesh(prob, data, ~, ub, Vb)
%cdata = coco_get_chart_data(chart, data.id);
stat='success';
pt=dde_point_from_x(ub,data.info.point,data.ipar);
ptv=dde_point_from_x(Vb,data.info.point,data.ipar);
ptnew=p_remesh(pt);
unew=dde_x_from_point(ptnew,data.ipar);
if any(unew~=ub)
    for i=length(ptv):-1:1
        ptvnew(i)=p_remesh(ptv(i),ptv(i).degree,ptnew.mesh);
    end
    Vnew=reshape(dde_x_from_point(ptvnew,data.ipar),size(Vb));
end
data.info.point=ptnew;
xtr=1:length(unew);
prob = coco_change_func(prob, data, 'u0', unew, 'vecs', Vnew); % Update data and solution
end
%%
function [data, res] = bddat(prob, data, command, varargin) %#ok<INUSL>
%BDDAT   store point struct in BD.

res = {};
switch command
  case 'init'
    res = coco_get_id(data.id,'point');
  case 'data'
    chart = varargin{1}; % Current chart
    res = dde_point_from_x(chart.x,data.info.point,data.ipar);
end
end
