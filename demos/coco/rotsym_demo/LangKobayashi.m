function y=LangKobayashi(E,Etau,n,p,ip)
%% r.h.s of Lang-Kobayashi equation
%
% $Id: LangKobayashi.m 20 2014-04-11 19:27:33Z jan.sieber $
%
Edot=(1+1i*p(ip.alpha))*n.*E+p(ip.eta)*exp(1i*p(ip.phi))*Etau;
ndot=p(ip.epsilon)*(p(ip.pump)-n-conj(E).*E.*(2*n+1));
y=cat(1,real(Edot),imag(Edot),real(ndot));
end
