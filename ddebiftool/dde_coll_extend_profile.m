function [dtext,err_ext]=dde_coll_extend_profile(dt,err_est)
%% function extending highest derivative by linear extrapolation  
% (used in auto_eqd)
%
% Input: dt: time increments, hd: error estimate
% Outputs: dtext,err_ext: extended time increments and error estimates
%
% $Id: dde_coll_extend_profile.m 308 2018-10-28 15:08:12Z jansieber $
%%
dtext=[2*dt(1)-dt(2),dt,2*dt(end)-dt(end-1)];
err_ext=[2*err_est(:,1)-err_est(:,2),err_est,2*err_est(:,end)-err_est(:,end-1)];
end