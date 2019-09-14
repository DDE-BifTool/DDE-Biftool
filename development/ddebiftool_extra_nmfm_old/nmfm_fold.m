function fold = nmfm_fold(funcs, fold, varargin)

sys_tau=funcs.sys_tau;
coord = fold.x;
par=fold.parameter;
tau=[0 par(sys_tau())];
r = length(tau); % Number of delays
xx = repmat(coord, 1, r); % All delayed vectors the same

Delta = ch_matrix(funcs,fold.x,par,0);
DDelta = ch_matrix(funcs,fold.x,par,0,'deri',1);

%% one dimensional case
if length(Delta) == 1  
    q=1;
    p=1;
else
 %% higher dimensional case
    q=fold.v;

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

%% critical normal forms coefficients
b=1/2*p'*funcs.sys_mfderi(xx,par, PHI, PHI);

fold.nmfm.b=b;
fold.nvec = q;


end