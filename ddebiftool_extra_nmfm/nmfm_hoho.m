function hoho = nmfm_hoho(funcs, point, varargin)
%% Compute normal form of Double Hopf point
%
% $Id: nmfm_hoho.m 314 2019-01-24 14:28:23Z mmbosschaert $
%
%%
default={'nullpoint',[],'print',0,'free_pars',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
kind = point.kind;
hoho = point;

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
[p1,q1]=nmfm_nullvector(funcs,point,lambda1,'nullpoint',options.nullpoint,pass_on{:});
[p2,q2]=nmfm_nullvector(funcs,point,lambda2,'nullpoint',[],pass_on{:});
%% abbreviate derivatives and characteristic matrix
F=nmfm_deriv_define(funcs,point,...
    'free_pars',options.free_pars,'print',options.print,pass_on{:});
Delta=@(z)ch_matrix(funcs,point.x,point.parameter,z);
par0=point.parameter(:)*0;
dev0=@(v)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))]);
devl=@(v,lambda)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',lambda);
devlt=@(v,lambda,t)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',lambda,'t',t);
%% Eigenfunctions
phi1 = devl(q1,lambda1);
phi1bar = nmfm_dev_conj(phi1);

phi2 = devl(q2,lambda2);
phi2bar = nmfm_dev_conj(phi2);

%% Normal form coefficients
% Quadratic center manifold
h1100 = dev0( Delta(0)\F.B(phi1,phi1bar));
h2000 = devl( Delta(2*lambda1)\F.B(phi1,phi1), 2*lambda1 );
h1010 = devl( Delta(lambda1+lambda2)\F.B(phi1,phi2), lambda1+lambda2 );
h1001 = devl( Delta(lambda1-lambda2)\F.B(phi1,phi2bar), lambda1-lambda2);
h1001bar = nmfm_dev_conj(h1001);
h0020 = devl( Delta(2*lambda2)\F.B(phi2,phi2), 2*lambda2);
h0011 = dev0( Delta(0)\F.B(phi2,phi2bar));

%% Cubic normal form
g2100 = (1/2)*p1*( 2*F.B(h1100,phi1) + F.B(h2000, phi1bar) + F.C(phi1,phi1,phi1bar) );
g1011 = p1*( F.B(h0011,phi1) + F.B(h1001,phi2) + F.B(h1010,phi2bar) + F.C(phi1,phi2,phi2bar) );
g1110 = p2*( F.B(h1001bar,phi1) + F.B(h1010,phi1bar)+ F.B(h1100,phi2) + F.C(phi1,phi1bar,phi2) );
g0021 = (1/2)*p2*( 2*F.B(h0011,phi2) + F.B(h0020,phi2bar) + F.C(phi2,phi2,phi2bar) );

theta0 = real(g1011)/real(g0021);
delta0 = real(g1110)/real(g2100);

%fprintf('theta(0) = %.10f, delta(0) = %.10f.\n', theta0, delta0);
%% Store normal form coefficients
hoho.nmfm.g2100 = g2100;
hoho.nmfm.g1011 = g1011;
hoho.nmfm.g1110 = g1110;
hoho.nmfm.g0021 = g0021;
hoho.nmfm.theta = theta0;
hoho.nmfm.delta = delta0;

hoho.nvec.p = [p1;p2];
hoho.nvec.q = [q1,q2];

%% parameter-related coefficients (Maikel, simplified by JS)
if isempty(options.free_pars)
    return
end
vp=eye(2);
for i=2:-1:1
    N(i)=dev0((Delta(0)\F.J1)*vp(:,i));
    M(1,i)=p1*(F.A1(phi1,vp(:,i))+F.B(phi1,N(i)));
    M(2,i)=p2*(F.A1(phi2,vp(:,i))+F.B(phi2,N(i)));
end
K=inv(real(M));
for i=2:-1:1
    h0000(i)=nmfm_dev_ax(-K(:,i),N);
end
for i=2:-1:1
    hoho.nmfm.b(i,1)=imag(p1*(F.A1(phi1,K(:,i))+F.B(phi1,h0000(i))));
    hoho.nmfm.b(i,2)=imag(p2*(F.A1(phi2,K(:,i))+F.B(phi2,h0000(i))));
end
hoho.nmfm.h0011=h0011;
hoho.nmfm.h0020=h0020;
hoho.nmfm.h2000=h2000;
hoho.nmfm.K=K;
hoho.nmfm.h0000=h0000;
end