function newpoint = nmfm_hopf(funcs, point, varargin)
%% Compute Lyapunov coefficient of Hopf point
%
% $Id: nmfm_hopf.m 315 2019-01-29 19:42:21Z jansieber $
%
%%
coord = point.x;
par = point.parameter;
kind = point.kind;
newpoint = point;

if ~strcmp(kind,'hopf')
    error('NMFM_HOPF: did not receive a hopf point, but %s, as argument.',kind);
end
if strcmp(point.flag,'zeho')
    newpoint.nmfm.L1=NaN;
    return
end
omega = abs(point.omega);
lambda0 = 1i*omega;
%% Compute nullvectors
[p,q,sg]=nmfm_nullvector(funcs,point,lambda0,varargin{:});
newpoint.nvec.omega = omega;
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
F=nmfm_deriv_define(funcs,point,'free_pars',[],varargin{:});
%% Construct additional characteristic matrices
Delta2 = ch_matrix(funcs, coord,par,2*lambda0);
Delta0 = ch_matrix(funcs, coord,par,0);
%% Implement normal form computations
par0=point.parameter(:)*0;
dev0=@(v)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))]);
devl=@(v,lambda)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',lambda);
phi = devl(q, lambda0);
phibar = nmfm_dev_conj(phi);

h20 = devl( Delta2\F.B(phi,phi), 2*lambda0);
h11 = dev0( 2*(Delta0\F.B(phi,phibar)));

c1 = (1/2)*p*( F.B(phibar,h20) + F.B(phi,h11) + F.C(phi,phi,phibar) );
L1 = real(c1)/(omega);
newpoint.nmfm.L1 = L1;
end
