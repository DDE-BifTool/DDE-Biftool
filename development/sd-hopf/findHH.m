function [ kappa_1,kappa_2,omega_1,omega_2, fv ] = findHH(j)
% findDH2(k1,k2,w1,w2) given approximate points (k1,k2)
% and corresponding w1 and w2 such that DDE has
% double Hopf point 
% (kappa_1,kappa_2) approx (k1,k2) with
% (omega_1,omega_2) approx (w1,w2)
% find
% (kappa_1,kappa_2) & (omega_1,omega_2)
%
% First three double Hopfs known to be close to
%
% (k1,k2)=(2.0805,3.7872) with (w1,w2) = (1.5820,2.4871) 
% (k1,k2)=(5.6086,2.6436) with (w1,w2) = (1.7710,6.6484)
% (k1,k2)=(9.2848,4.4039) with (w1,w2) = (1.9532,10.8246)

%% set other parameters
gam=4.75; epsi=1; a1=1.3; a2=6;


%% set (k1,k2,w1,w2) guess for required Hopf-Hopf point
if j==1
   k1=2.0805; k2=3.7872; w1=1.5856; w2=2.5363;
elseif j==2
   k1=5.6086; k2=2.6436; w1=1.7710; w2=6.6484;
elseif j==3
   k1=9.2848; k2=4.4039; w1=1.9532; w2=10.8246;
else
   disp(['HH',num2str(j),' not implemented'])
   stop
end

%% find the double Hopf point
options = optimset('TolFun',1e-14,'TolX',1e-14);
[x,fv]=fminsearch(@(x) dhres(x,a1,a2,gam,epsi),[k1,k2,w1,w2],options);
kappa_1=x(1);
kappa_2=x(2);
omega_1=max(x(3),x(4));   %  SWAP MAX AND MIN ON THESE TWO 
omega_2=min(x(3),x(4));   %  LINES TO SWAP ORDER OF THE OMEGAS

disp(['double Hopf Point found at kappa_1=',num2str(kappa_1,10),' kappa_2=',num2str(kappa_2,10),' omega_1=',num2str(omega_1,10),' omega_2=',num2str(omega_2,10)])
disp(['with error norm=',num2str(fv)])

end

function chr=chreal(k1,k2,w,a1,a2,gam)
   chr=gam+k1*cos(a1*w)+k2*cos(a2*w);
end

function chi=chimag(k1,k2,w,a1,a2,epsi)
   chi=epsi*w-k1*sin(a1*w)-k2*sin(a2*w);
end

function dh=dhres(x,a1,a2,gam,epsi)
  k1=x(1); k2=x(2); w1=x(3); w2=x(4);
  dh=norm([chreal(k1,k2,w1,a1,a2,gam); chimag(k1,k2,w1,a1,a2,epsi); chreal(k1,k2,w2,a1,a2,gam); chimag(k1,k2,w2,a1,a2,epsi)]);
end



