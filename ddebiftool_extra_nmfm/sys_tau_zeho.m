function tau=sys_tau_zeho(itau,xx,par,funcs)
%% delay. for Zero-Hopf bifurcation
%
% $Id$
%%
if ~funcs.tp_del
    tau=funcs.sys_tau();
else
n=size(xx,1)/4;
tau=funcs.sys_tau(itau,xx(1:n,:),par(1,1:end-1,:));
end
