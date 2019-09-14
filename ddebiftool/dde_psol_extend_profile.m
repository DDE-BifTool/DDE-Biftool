function [dtext,err_ext]=dde_psol_extend_profile(dt,err_est)
%% function extending highest derivative periodically
% (used in auto_eqd)
%
% Input: dt: time increments, hd: error estimate
% Outputs: dtext,err_ext: extended time increments and error estimates
%
% $Id: dde_psol_extend_profile.m 308 2018-10-28 15:08:12Z jansieber $
%%
dtext=[dt(end),dt,dt(1)];
err_ext=[err_est(:,end),err_est,err_est(:,1)];
end
