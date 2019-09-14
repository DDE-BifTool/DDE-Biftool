function [r,J]=sys_cond_fixperiod(point,ind_period)
%% fix period to given parameter (visible to user fcns)
% implemented as user condition
% $Id: sys_cond_fixperiod.m 369 2019-08-27 00:07:02Z jansieber $
%% 
r=point.parameter(ind_period)-point.period;
J=p_axpy(0,point,[]);
J.parameter(ind_period)=1;
J.period=-1;
end
