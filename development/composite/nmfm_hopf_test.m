function newpoint = nmfm_hopf(funcs,nfuncs, point, varargin)
%% Compute Lyapunov coefficient of Hopf point
%
% $Id: nmfm_hopf.m 185 2017-04-11 13:46:45Z mmbosschaert $
%
%%
coord = point.x;
par = point.parameter;
kind = point.kind;
n = length(coord);
newpoint = point;

if ~strcmp(kind,'hopf')
    error('NMFM_HOPF: did not receive a hopf point, but %s, as argument.',kind);
end
omega = abs(point.omega);
lambda0 = 1i*omega;
%% Compute nullvectors
[p,q,sg]=nmfm_nullvector(funcs,point,lambda0,varargin{:});
if sg
    newpoint.nvec.q = NaN(size(q));
    newpoint.nvec.p = NaN(size(p));
    newpoint.nmfm.L1= NaN;
    return
else
    newpoint.nvec.q = q;
    newpoint.nvec.p = p;
end
%% Select eigenvalue pair
if isempty(omega) || omega == 0
    warning('NMFM_HOPF:omega',...
        'NMFM_HOPF:  omega is empty or zero, returning L1 = NaN.');
    newpoint.nvec.q = NaN(size(q));
    newpoint.nvec.p = NaN(size(p));
    newpoint.nmfm.L1 = NaN;
    return
end
%% abbreviate derivative
B=@(v1,v2)nmfm_mfderiv(funcs,point,[v1,v2],varargin{:});
C=@(v1,v2,v3)nmfm_mfderiv(funcs,point,[v1,v2,v3],varargin{:});
Bn=@(v1,v2)nmfm_mfderiv_sav(nfuncs,point,[v1,v2],varargin{:});
Cn=@(v1,v2,v3)nmfm_mfderiv_sav(nfuncs,point,[v1,v2,v3],varargin{:});
%% Construct additional characteristic matrices
Delta2 = ch_matrix(funcs, coord,par,2*lambda0);
Delta0 = ch_matrix(funcs, coord,par,0);
%% Implement normal form computations
par0=point.parameter(:)*0;
phi = nmfm_dev_fun([q;par0],'lambda',lambda0);
phibar = nmfm_dev_conj(phi);

bnpp=Bn(phi,phi);
bpp=B(phi,phi);
bnppb=Bn(phi,phibar);
bppb=B(phi,phibar);
h20 = nmfm_dev_fun( [Delta2\B(phi,phi);par0], 'lambda',2*lambda0);
h11 = nmfm_dev_fun( [Delta0\B(phi,phibar);par0] );
h20n = nmfm_dev_fun( [Delta2\Bn(phi,phi);par0], 'lambda',2*lambda0);
h11n = nmfm_dev_fun( [Delta0\Bn(phi,phibar);par0] );

bph20=B(phibar,h20)
bph20n=Bn(phibar,h20n)
cppp=C(phi,phi,phibar);
cnppp=Cn(phi,phi,phibar);

c1 = (1/2)*p*( B(phibar,h20) + 2*B(phi,h11) + C(phi,phi,phibar) );
c1n = (1/2)*p*( Bn(phibar,h20) + 2*Bn(phi,h11) + Cn(phi,phi,phibar) );
L1 = real(c1)/(omega);
newpoint.nmfm.L1 = L1;
end
