function [hcli_br,hcli1,hcli2] = init_BT_Hom(funcs,bt,amplitude,TTolerance,finemsh_points)
%% Initialize branch for continuing the homoclinic curve
% emanating from the Bogdanov-Takens point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_BT_Hom.m 309 2018-10-28 19:02:42Z jansieber $
%%

if ~strcmp(bt.kind,'BT')
    display(bt.kind);
    error('init_BT_Hom: did not receive a bt point as argument.');
end

%% approximate homoclinic orbits
hcli1=hom_approx(funcs,bt,amplitude,TTolerance,finemsh_points);
hcli2=hom_approx(funcs,bt,amplitude*1.01,TTolerance,finemsh_points);
%% setup homoclinic branch
hcli_br=df_brnch(funcs,get_free_pars,'hcli');
hcli_br.point=hcli1;
hcli_br.point(2)=hcli2;

end

function hcli1 = hom_approx(funcs,bt,amplitude,TTolerance,finemsh_points)

sys_tau=funcs.sys_tau;

x=bt.x;
par=bt.parameter;
tau=[0 par(sys_tau())];
m=length(tau)-1;
xx=x(:,ones(m+1,1));

%% set coefficients
a=bt.nmfm.a;
b=bt.nmfm.b;
d=bt.nmfm.d;
e=bt.nmfm.e;
a2=bt.nmfm.a2;
b2=bt.nmfm.b2;

phi_0=bt.nmfm.phi_0;
phi_1=bt.nmfm.phi_1; % dependent on theta

h0001=bt.nmfm.h0001;
h1001=bt.nmfm.h1001;
h0002=bt.nmfm.h0002;
h2000=bt.nmfm.h2000; % dependent on theta
h0010=bt.nmfm.h0010; % dependent on theta

K10=bt.nmfm.K10;
K02=bt.nmfm.K02;
K01=bt.nmfm.K01;

%% saddle point approximation
eps=sqrt(amplitude*abs(a)/(6*sqrt(m+1)));
finemsh=0:1/(finemsh_points-1):1;

myarg = sqrt(TTolerance / (amplitude*sqrt(m+1)) );
T=1/eps*asech(myarg);
period = 2*T;

x0=xx(:,1);
theta=0;
x0 = x0 + eps^2 * ((10*b)/(7*a) * h0001 + 2/a * phi_0) ...
    +(phi_0*(-(2/7)*(5*a2*b+7*d)/a^2)...
    /a+(10/7)*h1001...
    *2*b/a^2 +(50/49)*h0002*b^2/a^2+(1/2)*h2000(theta)*...
    2^2/a^2+h0001*((100/49)*b2/a-(50/49)*a2...
    *b/a^2+(288/2401)*b^2/a^2- 4*e/(a*b) + (146/49)*d/a^2)*b/a-4*...
    h0010(theta)/a)*eps^4;

%% homoclinic orbit approximation
n=length(phi_0);
ups = zeros(n,length(finemsh));
for i=1:length(finemsh)
    t = (2*finemsh(i) - 1) * (T); % CONVERSION FROM [0,1] --> [-T,+T]
    
    theta=0;
    ups(:,i)=((10/7)*h0001*b/a + phi_0*(-6*sech(eps*t)^2+2)/a)*eps^2+...
        12*phi_1(theta)*sech(eps*t)^2*tanh(eps*t)*eps^3/a + (phi_0*(-(1/49)*(-210*...
        a2*b+18*b^2+ 147*d)*sech(eps*t)^2/a^2-(2/7)*(5*a2*b+7*d)/a^2)...
        /a- (72/7)*phi_1(theta)*b*sinh(eps*t)^2/(a^2*cosh(eps*t)^4)+(10/7)*h1001...
        *(-6*sech(eps*t)^2+2)*b/a^2 +(50/49)*h0002*b^2/a^2+(1/2)*h2000(theta)*...
        (-6*sech(eps*t)^2+2)^2/a^2+h0001*((100/49)*b2/a-(50/49)*a2...
        *b/a^2+(288/2401)*b^2/a^2- 4*e/(a*b) + (146/49)*d/a^2)*b/a-4*...
        h0010(theta)/a)*eps^4;
    
    ups(:,i) = ups(:,i) + bt.x;
end

%% parameter approximation
alpha=((10*b)/(7*a))*eps^2*K01+eps^4*((-4/a)*K10+((50*b^2)/(49*a^2))*K02+...
    b/a*K01*(1/a*(100/49*b2-4*e/4)+1/a^2*(-50/49*b*a2+288/2401*b^2+146/49*d)))+...
    bt.parameter(get_free_pars())';

%% setup homoclinic orbit
free_pars=get_free_pars();

hcli1.kind = 'hcli';
hcli1.parameter = bt.parameter;
hcli1.parameter(free_pars)=alpha;
hcli1.mesh = finemsh;
hcli1.degree = 3;
hcli1.profile = ups;
hcli1.period = period;
hcli1.x1=x0; % saddle
hcli1.x2=x0; % (homoclinic)
%
hcli1=p_tohcli(funcs,hcli1); % this adds the fields lambda_{v,w}, v,w, alpha,espilon

% correct homoclinic orbit
mhcli=df_mthod(funcs,'hcli');
[hcli1,s]=p_correc(funcs,hcli1,free_pars(1),[],mhcli.point); % correct
if ~s
    warning('Correction homoclinic orbit failed');
end

end
