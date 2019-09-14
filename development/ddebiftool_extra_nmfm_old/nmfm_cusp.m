function cusp = nmfm_cusp(funcs, cusp, varargin)

sys_tau=funcs.sys_tau;
coord = cusp.x;
par=cusp.parameter;
tau=[0 par(sys_tau())];
r = length(tau); % Number of delays
xx = repmat(coord, 1, r); % All delayed vectors the same

Delta = ch_matrix(funcs,cusp.x,par,0);
DDelta = ch_matrix(funcs,cusp.x,par,0,'deri',1);

%% one dimensional case
if length(Delta) == 1  
    q=1;
    p=1;
else
 %% higher dimensional case
    q=cusp.v;

    [V,D]=eig(Delta');
    [~,i2]=min(abs(diag(D)));
    p=V(:,i2);
end

%% normalize
alpha=p'*DDelta*q;
p=p/alpha;

%% eigenfunctions
phi = @(theta) q;
PHI = nmfm_handletomatrix(phi, -tau);

Bord = [Delta q; p' 0];
DinvD2f = Bord\[funcs.sys_mfderi(xx,par,PHI,PHI); 0];
DinvD2f= DinvD2f(1:2);

h2 = @(theta) -DinvD2f + p'*DDelta*(DinvD2f)*q;
H2 = nmfm_handletomatrix(h2, -tau);

%% critical normal forms coefficients
b=1/2*p'*funcs.sys_mfderi(xx,par, PHI, PHI);
c=1/6*p'*(3*funcs.sys_mfderi(xx,par, PHI, H2)+funcs.sys_mfderi(xx,par, PHI, PHI,PHI));

% b=b/alpha;
% c=c/alpha;

cusp.nmfm.b=b;
cusp.nmfm.c=c;
cusp.nvec = q;


end
