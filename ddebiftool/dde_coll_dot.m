function [y,J]=dde_coll_dot(point1,point2,varargin)
%% scalar product for coll points and its bifurcations
% $Id: dde_coll_dot.m 362 2019-07-14 15:49:40Z jansieber $
%%
default={'free_par_ind',1:length(point1.parameter),'period',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[p1,cv]=dde_coll_convert(point1);
p2=dde_coll_convert(point2);
ind=dde_ind_from_point(point1,options.free_par_ind);
[dum,W]=dde_coll_profile_dot(p1,p2,pass_on{:}); %#ok<ASGLU>
Jprof=W*p2.profile(:);
free_par_ind=[options.free_par_ind,cv.extra_ind.parameter];
p2.profile=reshape(Jprof,size(p1.profile));
x1=dde_x_from_point(p1,free_par_ind);
x2=dde_x_from_point(p2,free_par_ind);
if ~options.period
    x1(ind.period)=0;
    x2(ind.period)=0;
end
y=x1'*x2;
J=dde_point_from_x(x2,p_axpy(0,p1),free_par_ind);
J=dde_coll_convert(J,cv);
end
