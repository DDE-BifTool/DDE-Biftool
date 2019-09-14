function newpoint = nmfm_genh(funcs, point, varargin)
%% Compute 2nd Lyapunov coefficient for generalized Hopf point
%
% $Id: nmfm_genh.m 79 2015-01-02 18:42:50Z jan.sieber $
% 
%%
coord = point.x;
par = point.parameter;
kind = point.kind;
n = length(coord);
ii = sqrt(-1);
newpoint = point;

if ~strcmp(kind,'genh')
    display(kind);
    error('NMFM_GENH: did not receive a generalized hopf point as argument.');
end

%% Select eigenvalue pair
omega = point.omega;
if isempty(omega) || omega == 0
    fprintf('NMFM_GENH: omega is empty or zero, returning L1 = NaN.\n');
    newpoint.nmfm.L2 = NaN;
    return;
end

lambda0 = ii*omega;

%% Construct characteristic matrix
if ~funcs.tp_del
   taus = par(funcs.sys_tau());
   taus = [0, taus]; % First delay zero
   r = length(taus); % Number of delays
   xx = repmat(coord, 1, r); % All delayed vectors the same
else % state-dependent delays
   r = funcs.sys_ntau() + 1; % Number of delays
   xx = repmat(coord, 1, r); % All delayed vectors the same
   taus = zeros(1,r); % First delay zero
   for i = 2:r
      taus(i) = funcs.sys_tau(i-1,xx(:,1:i-1),par);
   end
end
Delta1 = ch_matrix(funcs,xx,par,lambda0);
Delta2 = ch_matrix(funcs,xx,par,2*lambda0);
Delta3 = ch_matrix(funcs,xx,par,3*lambda0);
Delta0 = ch_matrix(funcs,xx,par,0);
DDelta1 = ch_matrix(funcs,xx,par,lambda0,'deri',1);
DDelta2 = ch_matrix(funcs,xx,par,2*lambda0,'deri',1);
D2Delta1 = ch_matrix(funcs,xx,par,lambda0,'deri',2);

%% Compute nullvectors
if n == 1
    p0 = 1;
    q0 = 1;
else
    if nargin == 2 % use null vectors of point
        nullpoint = point;
    elseif nargin == 3 % use supplied null vectors
        nullpoint = varargin{1};
    end
    % Construct null vectors
    if isfield(nullpoint, 'nvec')
       if isfield(nullpoint.nvec,'p')
           pp0 = nullpoint.nvec.p; 
       else
           pp0 = [];
       end
       if isfield(nullpoint.nvec,'q') 
           qq0 = nullpoint.nvec.q; 
       else
           qq0 = []; 
       end
    else
       pp0 = []; qq0 = [];
    end
    [p0, q0] = nmfm_border(Delta1, pp0, qq0);
    p0 = p0/norm(p0);
    q0 = q0/norm(q0);
end

if isempty(p0) || isempty(q0)
   fprintf('NMFM_GENH: could not construct null vectors.\n');
   return;
end

%% Normalize eigenvectors
alpha = 1/sqrt(p0*DDelta1*q0);
p = alpha*p0;
q = alpha*q0;

%% Implement normal form computations
phi = @(theta) exp(lambda0*theta)*q;
phibar = @(theta) conj(exp(lambda0*theta)*q);

PHI = nmfm_handletomatrix(phi, -taus);
PHIBAR = nmfm_handletomatrix(phibar, -taus);

h20 = @(theta) exp(2*lambda0*theta)*(Delta2\funcs.sys_mfderi(xx,par,PHI,PHI));
H20 = nmfm_handletomatrix(h20,-taus);

h11 = @(theta) Delta0\funcs.sys_mfderi(xx,par, PHI, PHIBAR);
H11 = nmfm_handletomatrix(h11,-taus);

c1 = (1/2)*p*(funcs.sys_mfderi(xx,par, PHIBAR, H20) + ...
    2*funcs.sys_mfderi(xx,par, PHI, H11) + funcs.sys_mfderi(xx,par,PHI,PHI,PHIBAR));

h30 = @(theta) exp(3*lambda0*theta)*(Delta3\(3*funcs.sys_mfderi(xx,par,PHI,H20) +...
    funcs.sys_mfderi(xx,par,PHI,PHI,PHI)));
H30 = nmfm_handletomatrix(h30,-taus);

h21 = @(theta) exp(lambda0*theta)*(-2*c1)*((1/2)*p*D2Delta1*q-theta)*q;
H21 = nmfm_handletomatrix(h21, -taus);

h31 = @(theta) exp(2*lambda0*theta)*(Delta2\(funcs.sys_mfderi(xx,par,PHIBAR, H30) + ...
    3*funcs.sys_mfderi(xx,par,PHI,H21) + 3*funcs.sys_mfderi(xx,par,H20,H11) + ...
    3*funcs.sys_mfderi(xx,par,PHI,PHIBAR,H20) + ...
    3*funcs.sys_mfderi(xx,par,PHI,PHI,H11) + funcs.sys_mfderi(xx,par,PHI,PHI,PHI,PHIBAR))) - ...
    6*c1*(Delta2\((DDelta2-eye(n)-theta*Delta2)*h20(theta)));
H31 = nmfm_handletomatrix(h31,-taus);

h22 = @(theta) Delta0\(2*funcs.sys_mfderi(xx,par,PHIBAR,H21) + ...
    2*funcs.sys_mfderi(xx,par, H11, H11) + 2*funcs.sys_mfderi(xx,par,PHI,conj(H21)) +...
    funcs.sys_mfderi(xx,par,H20,conj(H20)) + funcs.sys_mfderi(xx,par,PHIBAR,PHIBAR,H20) +...
    funcs.sys_mfderi(xx,par,PHI,PHI,conj(H20)) + 4*funcs.sys_mfderi(xx,par,PHI,PHIBAR,H11) +...
    funcs.sys_mfderi(xx,par,PHI,PHI,PHIBAR,PHIBAR));
H22 = nmfm_handletomatrix(h22,-taus);

c2 = (1/12)*p*(6*funcs.sys_mfderi(xx,par,H11,H21) + 3*funcs.sys_mfderi(xx,par,conj(H21),H20) + ...
    funcs.sys_mfderi(xx,par,conj(H20),H30) + ...
    3*funcs.sys_mfderi(xx,par,PHI,H22) + 2*funcs.sys_mfderi(xx,par,PHIBAR,H31) + ...
    6*funcs.sys_mfderi(xx,par,PHIBAR,H20,H11) + 6*funcs.sys_mfderi(xx,par,PHI,H11,H11) + ...
    3*funcs.sys_mfderi(xx,par,PHI,H20,conj(H20)) + 6*funcs.sys_mfderi(xx,par,PHI,PHIBAR,H21) + ...
    3*funcs.sys_mfderi(xx,par,PHI,PHI,conj(H21)) + funcs.sys_mfderi(xx,par,PHIBAR,PHIBAR,H30) + ...
    6*funcs.sys_mfderi(xx,par,PHI,PHI,PHIBAR,H11) + ...
    3*funcs.sys_mfderi(xx,par,PHI,PHIBAR,PHIBAR,H20) + ...
    funcs.sys_mfderi(xx,par,PHI,PHI,PHI,conj(H20)) + ...
    funcs.sys_mfderi(xx,par,PHI,PHI,PHI,PHIBAR,PHIBAR));


L2 = real(c2)/(omega);

%% Attach result to point structure
if ~isfield(newpoint.nmfm,'L1') || isempty(newpoint.nmfm.L1) || isnan(newpoint.nmfm.L1)
    newpoint.nmfm.L1 = real(c1)/omega;
end

newpoint.nmfm.L2 = L2;
newpoint.nvec.p = p;
newpoint.nvec.q = q;

end
