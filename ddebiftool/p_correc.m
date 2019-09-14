function [point,success]=p_correc(funcs,point,free_par,step_cnd,method,remesh_flag,...
    previous,varargin)
%% correct point using Newton iteration
% function [point,success]=p_correc(point0,free_par,step_cnd,method,adapt,
%                                 previous)
% INPUT:
%   funcs problem functions
%   point0 initial point guess
%   free_par free parameter numbers in N^d
%       step_cnd steplength condition(s) as point(s)
%   method method parameters
%   remeshflag if zero or absent, do not adapt mesh; if one, always adapt
%       previous (optional) previously computed branch point (used, in case
%            of periodic solutions or connecting orbits, to
%            minimize phase shift)
% OUTPUT:
%       point final solution
%       success nonzero for successfull computation
%
% (c) DDE-BIFTOOL v. 2.02, 16/6/2002
%
% $Id: p_correc.m 348 2019-06-19 13:09:01Z jansieber $
%
%%
%% optional parameters:
if nargin<=5
    remesh_flag=0;
end
if nargin<=6
    previous=[];
end
[f,x0,data]=p_correc_setup(funcs,point,free_par,method,...
    'remesh_flag',remesh_flag,'previous',previous,'step_cnd',step_cnd,varargin{:});
%% perform Newton-Raphson iterations:
[x,success]=dde_nsolve(f,x0,data.method);
data.point=dde_point_from_x(x,data.point,data.free_par);
point=dde_postprocess(data.point,data);
%% recorrect with adapted mesh if necessary
if success && isfield(method,'remesh') && method.remesh
    if isfield (point,'mesh') && ~isempty(point.mesh)
        ma=method.adapt_mesh_after_correct;
        if remesh_flag==1 || ( remesh_flag>1 && mod(remesh_flag,ma)==0 )
            % do not adapt mesh when p_nr=0
            % adapt mesh when p_nr=1
            % adapt mesh when p_nr>1 & mod(p_nr,ma)=0
            if success % adapt & correct again
                method2=method;
                method2.adapt_mesh_before_correct=1;
                method2.adapt_mesh_after_correct=0;
                [point,success]=p_correc(funcs,point,free_par,step_cnd,method2,2,previous,...
                    varargin{:});
            end
        end
    end
end

end
