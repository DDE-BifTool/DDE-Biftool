function [trfuncs,trbranch,psol1] = init_HH_NS1(funcs,HH,eps,varargin)
%% Initialize branch for continuing the first Neimark-Sacker curve
% emanating from the the Hopf-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_HH_NS1.m 309 2018-10-28 19:02:42Z jansieber $
%%

%% approximate first two cycles
psol1=cycle_approx(HH,eps);
psol2=cycle_approx(HH,eps*1.01);

%% Initialize Neimark-Sacker branch
[trfuncs,trbranch]=init_NS_br(funcs,psol1,psol2,varargin{:});
end

function psol = cycle_approx(HH,eps)
%% normal form coefficients
g1110=HH.nmfm.g1110;
g2100=HH.nmfm.g2100;
H0011=HH.nmfm.H0011;
H2000=HH.nmfm.H2000;

b11=HH.nmfm.b11;
b12=HH.nmfm.b12;

K10=HH.nmfm.K10;
K01=HH.nmfm.K01;
H000010=HH.nmfm.H000010;
H000001=HH.nmfm.H000001;

q1=HH.nvec.q1;
%% approxiamtion to the cycle
newbase = HH.x + (H0011(:,1) -real(g1110)*H000001(:,1)-real(g2100)*H000010(:,1))*eps^2;
n = length(HH.x);
x0 = zeros(n,81);
for i=1:81
    zz = -exp(sqrt(-1.0)*(2*pi*(i-1)/80+pi/2));
    x0(:,i) = newbase + real(2*q1*zz*eps + H2000(:,1)*zz^2*eps^2);
end
%% approximation to the period
beta1=-real(g2100)*eps^2;
beta2=-real(g1110)*eps^2;
T1=2*pi/(HH.omega1+b11*beta1+b12*beta2+imag(g2100)*eps^2);

%% setup periodic orbit
psol.kind='psol';
psol.parameter=HH.parameter;
pm=-(real(g2100)*K10+real(g1110)*K01)*eps^2;
psol.parameter(get_free_pars())=...
    real(psol.parameter(get_free_pars())+pm');
psol.mesh=linspace(0,1,81);
psol.degree=4;
psol.profile=real(x0);
psol.period=real(T1);
end
