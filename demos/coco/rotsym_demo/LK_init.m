function [E0,n0,par]=LK_init(par,ip)
%% obtain initial value for rotating wave for given alpha, pump, eta and tau
%
% $Id: LK_init.m 20 2014-04-11 19:27:33Z jan.sieber $
%
n0=par(ip.eta)/2;
omega0=n0*par(ip.alpha)-sqrt(par(ip.eta)^2-n0^2);
phi0=angle(n0*(1+1i*par(ip.alpha))-1i*omega0)+pi+omega0*par(ip.tau);
phi0=mod(phi0,2*pi);
E0=sqrt((par(ip.pump)-n0)/(2*n0+1));
E0=[E0;0];
par(ip.phi)=phi0;
par(ip.omega)=omega0;
end
