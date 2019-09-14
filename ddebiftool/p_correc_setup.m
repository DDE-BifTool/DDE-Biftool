function [f,x0,data]=p_correc_setup(funcs,point,free_par,method,varargin)
%% set up f and initial guess x0,p reprocess and remesh point if required
%
% INPUT:
%
%  * funcs problem functions
%  * point initial point guess
%  * free_par free parameter numbers in N^d
%  * method: point method parameters
%  * step_cnd (default  empty) steplength condition(s) as point(s)
%  * remesh_flag (default 0) if zero or absent, do not adapt mesh; if one, always adapt
%  * previous (default empty) previously computed point along branch (for
%  example for periodic solutions or connecting orbits, to minimize phase
%  shift)
%
% OUTPUT:
%  * f function of type [r,J]=f(x) returning residual
%  * x0: vector, initial guess corresponding to point
%  * data: structure with remaining info (preprocessed free_par, method,
%  etc)
% 
% $Id: p_correc_setup.m 352 2019-06-19 16:03:03Z jansieber $
%%
default={'remesh_flag',false,'step_cnd',repmat(point,0,1),...
    'previous',repmat(point,0,1),'extracolumns',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% replace old stability and other info by empty fields if present:
point=dde_trim_point(feval(['dde_',point.kind,'_create'],point),point);

%% preprocess point, stepcond and functions
data=struct('funcs',funcs,'free_par',free_par,'method',method,...
    'step_cnd',options.step_cnd,'point',point,...
    'previous',options.previous,'extracolumns',options.extracolumns);
if isfield(method,'preprocess') && ~isempty(method.preprocess)
    data=feval(method.preprocess,data);
end
%% remesh if required
if isfield(data.method,'remesh') && ...
        data.method.remesh && ...
        options.remesh_flag>1 && ...
        mod(options.remesh_flag,data.method.adapt_mesh_before_correct)==0
    data.point=p_remesh(data.point);
    if ~isempty(data.previous)
        data.previous=dde_trim_point(p_remesh(data.previous,data.point),point);
    end
    for i=1:length(data.step_cnd)
        data.step_cnd(i)=dde_trim_point(p_remesh(data.step_cnd(i),data.point),point);
    end
end
%% prepare f and x0 for Newton-Raphson iterations:
f=@(x)p_correc_rhs(data.funcs,data.method,data.point,data.free_par,...
    'x',x,'pref',data.previous,'step_cond',data.step_cnd,...
    'extrapar',data.extracolumns,pass_on{:});
if nargout>1
    x0=dde_x_from_point(data.point,data.free_par);
end
end
