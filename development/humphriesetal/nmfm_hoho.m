function newpoint = nmfm_hoho(funcs, point, varargin)
%% Compute normal form of Double Hopf point
%
% $Id: nmfm_hoho.m 152 2017-02-17 22:57:53Z jansieber $
%
%%
default={'nullpoint',[],'debug',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
kind = point.kind;
newpoint = point;

if ~strcmp(kind,'hoho')
   display(kind);
   error('NMFM_HOHO: did not receive a hopf-hopf point as argument.');
end

%% Select eigenvalue pairs
omega1 = point.nvec.omega(1);
omega2 = point.nvec.omega(2);
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
if options.debug
    Delta=@(z,k)ch_matrix(funcs,point.x,point.parameter,z,'deri',k);
    p1q1=p1*Delta(lambda1,1)*q1;
    p2q2=p2*Delta(lambda2,1)*q2;
    small=@(x)norm(x)<1e-12;
    assert(small(p1q1-1));
    assert(small(p2q2-1));
end
%% abbreviate derivatives and characteristic matrix
B=@(v1,v2)nmfm_mfderiv(funcs,point,[v1,v2],pass_on{:});
C=@(v1,v2,v3)nmfm_mfderiv(funcs,point,[v1,v2,v3],pass_on{:});
Delta=@(z)ch_matrix(funcs,point.x,point.parameter,z);
%% Eigenfunctions
phi1 = nmfm_dev_fun(q1,'lambda',lambda1);
phi1bar = nmfm_dev_conj(phi1);

phi2 = nmfm_dev_fun(q2,'lambda',lambda2);
phi2bar = nmfm_dev_conj(phi2);

%% Normal form coefficients
% Quadratic center manifold
h1100 = nmfm_dev_fun( Delta(0)\B(phi1,phi1bar) );
h2000 = nmfm_dev_fun( Delta(2*lambda1)\B(phi1,phi1), 'lambda',2*lambda1 );
h1010 = nmfm_dev_fun( Delta(lambda1+lambda2)\B(phi1,phi2), 'lambda',lambda1+lambda2 );
h1001 = nmfm_dev_fun( Delta(lambda1-lambda2)\B(phi1,phi2bar), 'lambda',lambda1-lambda2);
h1001bar = nmfm_dev_conj(h1001);
h0020 = nmfm_dev_fun( Delta(2*lambda2)\B(phi2,phi2), 'lambda',2*lambda2);
h0011 = nmfm_dev_fun( Delta(0)\B(phi2,phi2bar) );

%% Cubic normal form
g2100 = (1/2)*p1*( 2*B(h1100,phi1) + B(h2000, phi1bar) + C(phi1,phi1,phi1bar) );
g1011 = p1*( B(h0011,phi1) + B(h1001,phi2) + B(h1010,phi2bar) + C(phi1,phi2,phi2bar) );
g1110 = p2*( B(h1001bar,phi1) + B(h1010,phi1bar)+B(h1100,phi2) + C(phi1,phi1bar,phi2) );
g0021 = (1/2)*p2*( 2*B(h0011,phi2) + B(h0020,phi2bar) + C(phi2,phi2,phi2bar) );

theta0 = real(g1011)/real(g0021);
delta0 = real(g1110)/real(g2100);

%fprintf('theta(0) = %.10f, delta(0) = %.10f.\n', theta0, delta0);

%% Store normal form coefficients
newpoint.nmfm.g2100 = g2100;
newpoint.nmfm.g1011 = g1011;
newpoint.nmfm.g1110 = g1110;
newpoint.nmfm.g0021 = g0021;
newpoint.nmfm.theta = theta0;
newpoint.nmfm.delta = delta0;

newpoint.nvec.p = [p1,p2];
newpoint.nvec.q = [q1,q2];

end