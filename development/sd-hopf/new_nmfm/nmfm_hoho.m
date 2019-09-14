function newpoint = nmfm_hoho(funcs, point, varargin)
%% Compute normal form of Double Hopf point
%
% $Id: nmfm_hoho.m 79 2015-01-02 18:42:50Z jan.sieber $
%
%%
default={'nullpoint',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
coord = point.x;
par = point.parameter;
kind = point.kind;
n = length(coord);
newpoint = point;

if ~strcmp(kind,'hoho')
   display(kind);
   error('NMFM_HOHO: did not receive a hopf-hopf point as argument.');
end

%% Select eigenvalue pairs
omega1 = point.omega1;
omega2 = point.omega2;
if isempty(omega1) || omega1 == 0
   fprintf('NMFM_HOHO: omega1 is empty or zero, unable to compute normal form.\n');
   return;
end
if isempty(omega2) || omega2 == 0
   fprintf('NMFM_HOHO: omega2 is empty or zero, unable to compute normal form.\n');
   return;
end

lambda1 = 1i*omega1;
lambda2 = 1i*omega2;
%% Compute nullvectors
[p1,q1]=nmfm_nullvector(funcs,point,lambda1,'nullpoint',options.nullpoint);
[p2,q2]=nmfm_nullvector(funcs,point,lambda2);
%% abbreviate derivative
mfd=@(dev)nmfm_mfderiv(funcs,point,dev,pass_on{:});
Delta=@(z)ch_matrix(funcs,coord,par,z);
%% Eigenfunctions
phi1 = nmfm_dev_fun(q1,'lambda',lambda1);
phi1bar = nmfm_dev_conj(phi1);

phi2 = nmfm_dev_fun(q2,'lambda',lambda2);
phi2bar = nmfm_dev_conj(phi2);

% Normal form coefficients
% Quadratic center manifold
h1100base=Delta(0)\mfd([phi1,phi1bar]);
h1100 = nmfm_dev_fun(h1100base);
h2000base=Delta(2*lambda1)\mfd([phi1,phi1]);
h2000 = nmfm_dev_fun(h2000base,'lambda',2*lambda1);
h1010base=Delta(lambda1+lambda2)\mfd([phi1,phi2]);
h1010 = nmfm_dev_fun(h1010base,'lambda',lambda1+lambda2);
h1001base=Delta(lambda1-lambda2)\mfd([phi1,phi2bar]);
h1001 = nmfm_dev_fun(h1001base,'lambda',lambda1-lambda2);
h0020base=Delta(2*lambda2)\mfd([phi2,phi2]);
h0020 = nmfm_dev_fun(h0020base,'lambda',2*lambda2);
h0011base=Delta(0)\mfd([phi2,phi2bar]);
h0011 = nmfm_dev_fun(h0011base);

% Cubic normal form
g2100 = (1/2)*p1*(2*mfd([h1100,phi1]) + mfd([h2000, phi1bar]) + mfd([phi1,phi1,phi1bar]));
g1011 = p1*(mfd([h0011,phi1])+mfd([h1001,phi2])+mfd([h1010,phi2bar])+mfd([phi1,phi2,phi2bar]));
g1110 = p2*(mfd([nmfm_dev_conj(h1001),phi1]) + mfd([h1010,phi1bar])+ ...
   mfd([h1100,phi2]) + mfd([phi1,phi1bar,phi2]));
g0021 = (1/2)*p2*(2*mfd([h0011,phi2]) + mfd([h0020,phi2bar]) + mfd([phi2,phi2,phi2bar]));

theta0 = real(g1011)/real(g0021);
delta0 = real(g1110)/real(g2100);

%fprintf('theta(0) = %.10f, delta(0) = %.10f.\n', theta0, delta0);

% Store normal form coefficients
newpoint.nmfm.g2100 = g2100;
newpoint.nmfm.g1011 = g1011;
newpoint.nmfm.g1110 = g1110;
newpoint.nmfm.g0021 = g0021;
newpoint.nmfm.theta = theta0;
newpoint.nmfm.delta = delta0;

newpoint.nvec.p = p1;
newpoint.nvec.q = q1;

end