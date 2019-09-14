function [ns_pfuncs,ns_br,psol1] = init_HT_NS1(funcs,ht,eps,varargin)
%% Initialize branch for continuing the first Neimark-Sacker curve
% emanating from the Hopf-Transcritical point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_HT_NS1.m 314 2019-01-24 14:28:23Z mmbosschaert $
%%

%% approximate first two cycles
psol1=cycle_approx(ht,eps);
psol2=cycle_approx(ht,eps*1.01);

%% Initialize Neimark-Sacker branch
[ns_pfuncs,ns_br] = init_NS_br(funcs,psol1,psol2,varargin{:});

end

function psol = cycle_approx(ht,eps)
%% normal form coefficients
g021=ht.nmfm.g021; g200=ht.nmfm.g200; g110=ht.nmfm.g110;
g011=ht.nmfm.g011; g210=ht.nmfm.g210; g300=ht.nmfm.g300;
g111=ht.nmfm.g111;

% K10=ht.nmfm.K10; K01=ht.nmfm.K01;
K10=ht.nmfm.K(:,1); K01=ht.nmfm.K(:,2);

% H20000=ht.nmfm.H200(:,1);
% H01100=ht.nmfm.H011(:,1);
% H02000=ht.nmfm.H020(:,1);
% H11000=ht.nmfm.H110(:,1);
% H10010=ht.nmfm.H10010(:,1);
% H10001=ht.nmfm.H10001(:,1);
% H01010=ht.nmfm.H01010(:,1);
% H01001=ht.nmfm.H01001(:,1);

n=length(ht.x);
get_n=@(x)x(1:n);
dev_eval=@(x)get_n(nmfm_dev_call(x,0));

H10010=dev_eval(ht.nmfm.h100mu(1));
H10001=dev_eval(ht.nmfm.h100mu(2));
H01010=dev_eval(ht.nmfm.h010mu(1)); % safely real
H01001=dev_eval(ht.nmfm.h010mu(2)); % safely real
H01100=real(dev_eval(ht.nmfm.h011)); % should be real anyway, but may be affected by roundoff
H02000=dev_eval(ht.nmfm.h020);
H20000=dev_eval(ht.nmfm.h200);
H11000=dev_eval(ht.nmfm.h200);

q0=ht.nvec.q0;
q1=ht.nvec.q1;

omega1=ht.nmfm.omega1;
omega2=ht.nmfm.omega2;

%% parmeters
beta1=2*sqrt(g011*g200).*eps-(g300*g011/g200+g111)*eps.^2;
beta2=real(g110)*sqrt(g011/g200).*eps+...
    (-g011*real(g210)/g200-g300*g011*real(g110)/(2*g200^2)...
    + real(g110)*real(g021)/g200-real(g021))*eps.^2;

%% period
z0=-sqrt(g011/g200)*eps-(2*g200*real(g021)-g300*g011)/(2*g200^2)*eps.^2;
T=2*pi/(ht.omega+omega1*beta1+omega2*beta2+...
    imag(g110)*z0+imag(g021)*eps^2);

%% cycle
nphase=1;
n = length(ht.x);
x0 = zeros(n,80*nphase+1);
for i=1:81
    psi=2*pi*(i-1)/80;
    x0(:,i) = ht.x + (-sqrt(g011/g200)*q0+2*real(exp(1i*psi)*q1))*eps ...
        + ( -(2*real(g021)*g200-g300*g011)/(2*g200^2)*q0 ...
        +2*real(g110)*sqrt(g011/g200)*real(exp(1i*psi)*H01001) ...
        + real(exp(2*1i*psi)*H02000) + H01100  ...
        - 2*sqrt(g011/g200)*real(exp(1i*psi)*H11000)...
        - real(g110)*g011/g200*H10001 ...
        + 4*sqrt(g011*g200)*real(exp(1i*psi)*H01010)...
        - 2*g011*H10010 + g011/(2*g200)*H20000)*eps^2;
end

%% setup periodic orbit
psol.kind='psol';
psol.parameter=ht.parameter;
psol.parameter(get_free_pars)=real(ht.parameter(get_free_pars)...
    +beta1*K10'+beta2*K01')';
psol.mesh=linspace(0,1,80*nphase+1);
psol.degree=4;
psol.profile=real(x0);
psol.period=real(T);
end
