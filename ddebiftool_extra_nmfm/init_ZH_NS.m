function [ns_pfuncs,ns_br,psol1] = init_ZH_NS(funcs,ZH,eps,varargin)
%% Initialize branch for continuing the Neimark-Sacker curve
% emanating from the Fold-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_ZH_NS.m 314 2019-01-24 14:28:23Z mmbosschaert $
%%

%% check condition for Neimark-Sacker bifurcation
if ZH.nmfm.g011*real(ZH.nmfm.g110)>0
    error('Condition for Neimark-Sacker bifurcation not met!')
end

%% approximate first two cycles
psol1=cycle_approx(ZH,eps);
psol2=cycle_approx(ZH,eps*1.01);

%% Initialize Neimark-Sacker branch
[ns_pfuncs,ns_br] = init_NS_br(funcs,psol1,psol2,varargin{:});

end

function psol = cycle_approx(ZH,eps,varargin)
% normal form coefficients
g111=ZH.nmfm.g111;
g021=ZH.nmfm.g021;
g200=ZH.nmfm.g200;
g110=ZH.nmfm.g110;
g011=ZH.nmfm.g011;

K10=ZH.nmfm.K(:,1);
K01=ZH.nmfm.K(:,2);

H00010=ZH.nmfm.h000mu(:,1);
H00001=ZH.nmfm.h000mu(:,2);
H011=ZH.nmfm.h011;
H020=ZH.nmfm.h020;

q0=ZH.nvec.q0;
q1=ZH.nvec.q1;

omega1=ZH.nmfm.omega1;
omega2=ZH.nmfm.omega2;

%% parameters
beta1=-g011*eps^2;
beta2 =(2*(real(g110)-g200)*real(g021)+real(g110)*g111)/(2*g200)*eps^2;
KK=ZH.parameter(get_free_pars())'+(beta1*K10+beta2*K01);

%% period
z0=-(2*real(g021)+g111)/(2*g200)*eps^2;
T=2*pi/(ZH.omega+omega1*beta1+omega2*beta2+imag(g021)*eps^2+imag(g110)*z0);

%% Limit cycle
n = length(ZH.x);
newbase = ZH.x+(-g011*H00010.v(1:n,1)*eps^2+...
    beta2*H00001.v(1:n,1)+z0*q0 + H011.v(1:n,1)*eps^2);
nphase=1;
x0 = zeros(n,80*nphase+1);
for i=1:81
    zz = -exp(sqrt(-1.0)*(2*pi*(i-1)/80+pi/2));
    x0(:,i) = newbase + real(2*q1*zz*eps + H020.v(1:n,1)*zz^2*eps^2);
end

%% setup periodic orbit
psol.kind='psol';
psol.parameter=ZH.parameter;
psol.parameter(get_free_pars)=real(KK)';
psol.mesh=linspace(0,1,80*nphase+1);
psol.degree=4;
psol.profile=real(x0);
psol.period=real(T);
end
