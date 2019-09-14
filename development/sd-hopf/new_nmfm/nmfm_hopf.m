function newpoint = nmfm_hopf(funcs, point, varargin)
%% Compute Lyapunov coefficient of Hopf point
%
% $Id: nmfm_hopf.m 79 2015-01-02 18:42:50Z jan.sieber $
%
%%
default={'nullpoint',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
coord = point.x;
par = point.parameter;
kind = point.kind;
n = length(coord);
newpoint = point;

if ~strcmp(kind,'hopf')
    error('NMFM_HOPF: did not receive a hopf point, but %s, as argument.',kind);
end

%% Select eigenvalue pair
omega = abs(point.omega);
if isempty(omega) || omega == 0
    warning('NMFM_HOPF:omega',...
        'NMFM_HOPF:  omega is empty or zero, returning L1 = NaN.');
    newpoint.nmfm.L1 = NaN;
    return
end

lambda0 = 1i*omega;
%% Compute nullvectors
[p,q]=nmfm_nullvector(funcs,point,lambda0,'nullpoint',options.nullpoint);
%% Construct additional characteristic matrices
Delta2 = ch_matrix(funcs, coord,par,2*lambda0);
Delta0 = ch_matrix(funcs, coord,par,0);

%% abbreviate derivative
mfderiv=@(dev)nmfm_mfderiv(funcs,point,dev,pass_on{:});
%% Implement normal form computations
phi = nmfm_dev_fun(q,'lambda',lambda0);
phibar = nmfm_dev_conj(phi);

D2_by_B_PHI_PHI=Delta2\mfderiv([phi,phi]);
D0_by_B_PHI_PHIbar=Delta0\mfderiv([phi,phibar]);

h20 = nmfm_dev_fun(D2_by_B_PHI_PHI,'lambda',2*lambda0);
h11 = nmfm_dev_fun(D0_by_B_PHI_PHIbar);

c1 = (1/2)*p*(mfderiv([phibar,h20])+...
    2*mfderiv([phi,h11]) +...
    mfderiv([phi,phi,phibar]));

L1 = real(c1)/(omega);

newpoint.nmfm.L1 = L1;
newpoint.nvec.p = p;
newpoint.nvec.q = q;

end
