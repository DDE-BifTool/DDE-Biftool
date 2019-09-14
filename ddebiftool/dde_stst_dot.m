function [y,J]=dde_stst_dot(point1,point2,varargin)
%% scalar product for stst points and its bifurcations
% $Id: dde_stst_dot.m 348 2019-06-19 13:09:01Z jansieber $
%%
default={'free_par_ind',1:length(point1.parameter)};
options=dde_set_options(default,varargin);
x1=dde_x_from_point(point1,options.free_par_ind);
x2=dde_x_from_point(point2,options.free_par_ind);
y=x1'*x2;
template=p_axpy(0,point2);
J=dde_point_from_x(x2,template,options.free_par_ind);
end

