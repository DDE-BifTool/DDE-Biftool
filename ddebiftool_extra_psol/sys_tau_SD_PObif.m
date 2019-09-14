function tau=sys_tau_SD_PObif(itau,x,p,ip,funcs)
%% delays of extended systems for fold, torus & period doubling
%
% $Id: sys_tau_SD_PObif.m 369 2019-08-27 00:07:02Z jansieber $
%
itau=mod(itau-1,ip.orig_ntau)+1;
tau=funcs.sys_tau(itau,x(1:ip.dim,1:itau,:),p(1,1:ip.nuserpar));
end
