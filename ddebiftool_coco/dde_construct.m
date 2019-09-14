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
%$Id: dde_construct.m 362 2019-07-14 15:49:40Z jansieber $
%% optional inputs
default={...
    'id','',...
    'funcs',[],...
    'info',[],...
    'pnames',{},...
    'parbasename','par',...
    'tangent',[],...
    'extpars',[],...
    'bd_point',true,...
    'fdf',true,'monitors',...
    {'xnorm',@xnorm_monitor}};
[options,pass_on]=dde_set_options(default,varargin,'pass_on','data');
%% construction of initial solution
% possibly redefine the funcs structure, create initial solution and
% parameters stored in typical "branch" structure of DDE-Biftool branch is
% supposed to contain one point and (optional) one tangent. This function
% calls, if present, dde_kind2type, where type is input, kind is field in
% options.info.point (or options.info). If not present it calls SetupType
% (Type capitalized) from ddebiftool_utilities.
options=add_extpars(options); % flag internal parameters as free if required
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
if options.fdf
    fdf={'f+df'};
else
    fdf={};
end
id=options.id;
data.id=id;
data.userfuncs=options.funcs;
data.monitors=options.monitors;
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
data.extpars=unique([data.info.parameter.free,options.extpars]);
data.userpars=setdiff(data.ipar,data.extpars);
if ~isempty(data.extpars)
    data.info.method.point.extra_condition=true;
end
data=coco_func_data(data);
prob=coco_add_func(prob,fid,@dde_rhs,data,'zero','u0',u0,t_info{:}, 'remesh', @dde_remesh,fdf{:});
uidx=coco_get_func_data(prob,fid,'uidx');
mon_names=coco_get_id(id,data.monitors(:,1));
monid=coco_get_id(id,'monitor');
prob=coco_add_func(prob,monid,@dde_monitors,struct('id',fid),'regular',mon_names,'uidx',uidx);
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
    if isfield(info,'kind') && strcmp(info.kind,'stst')
        info_ind={'point',info};
    else
        info_ind={'point',info.point};
    end
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
function [data, y] = dde_monitors(prob, data, u)
ptdata=coco_get_func_data(prob,data.id,'data');
pt=dde_point_from_x(u,ptdata.info.point(1),ptdata.ipar);
for i=size(ptdata.monitors,1):-1:1
    y(i)=ptdata.monitors{i,2}(pt);
end
end
%%
function y=xnorm_monitor(pt)
if ~isempty(ptdata.monitors)
    y=ptdata.monitors(pt);
elseif isfield(pt,'x')
    y=norm(pt.x,2);
elseif isfield(pt,'profile')
    y=sqrt(p_dot(pt,pt));
end
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
%%
function pnames=fill_pnames(names,base,ind)
pnames=arrayfun(@(i)[base,num2str(i)],ind,'uniformoutput',false);
if isstruct(names)
    fnames=fieldnames(names);
    cind=struct2cell(names);
    if ~isempty(ind)
        [~,ism]=ismember([cind{:}],ind);
        pnames(ism)=fnames;
    else
        pnames=fnames;
    end
else
    pnames(1:length(names))=names;
end
%pnames=coco_get_id(parid,pnames);
end
%%
function options=add_extpars(options)
%% are internal parameters free?
if  ischar(options.extpars)
    options.extpars={options.extpars};
end
if iscell(options.extpars)
    pnames=fill_pnames(options.pnames,options.parbasename,[]);
    extparnames=options.extpars;
    options.extpars=NaN(1,length(options.extpars));
    for i=length(options.extpars):-1:1
        options.extpars(i)=find(strcmp(extparnames{i},pnames));
    end
end
if isfield(options.info,'method') %info is a branch
    options.info.parameter.free=options.extpars;
end
end