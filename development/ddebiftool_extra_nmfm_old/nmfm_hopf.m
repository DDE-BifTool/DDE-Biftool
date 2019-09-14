function newpoint = nmfm_hopf(funcs, point, varargin)
%% Compute Lyapunov coefficient of Hopf point
%
% $Id: nmfm_hopf.m 151 2017-02-14 11:38:46Z jansieber $
%
%%

coord = point.x;
par = point.parameter;
kind = point.kind;
n = length(coord);
ii = sqrt(-1);
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

Delta1 = ch_matrix(funcs, xx,par,lambda0);
Delta2 = ch_matrix(funcs, xx,par,2*lambda0);
Delta0 = ch_matrix(funcs, xx,par,0);
DDelta1 = ch_matrix(funcs,xx,par,lambda0,'deri',1);

%% Compute nullvectors
if n == 1
    p0 = 1;
    q0 = 1;
else
    if nargin == 2 % use null vectors of point
        nullpoint = point;
    elseif nargin == 3 % use supplied null vectors
        nullpoint = varargin{1};
    elseif nargin == 4
        nullpoint = varargin{2};
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
   warning('NMFM_HOPF:nullvectors',...
       'NMFM_HOPF: null vectors are empty, returning L1 = NaN.\n');
   newpoint.nmfm.L1 = NaN;
   return
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

D2_by_B_PHI_PHI=Delta2\funcs.sys_mfderi(xx,par,PHI,PHI);
D0_by_B_PHI_PHIbar=Delta0\funcs.sys_mfderi(xx,par,PHI,PHIBAR);

h20 = @(theta) exp(2*lambda0*theta)*D2_by_B_PHI_PHI;
h11 = @(theta) D0_by_B_PHI_PHIbar;
H20 = nmfm_handletomatrix(h20,-taus);
H11 = nmfm_handletomatrix(h11,-taus);

c1 = (1/2)*p*(funcs.sys_mfderi(xx,par, PHIBAR, H20) + ...
    2*funcs.sys_mfderi(xx,par, PHI, H11) +...
    funcs.sys_mfderi(xx,par,PHI,PHI,PHIBAR));

L1 = real(c1)/(omega);

newpoint.nmfm.L1 = L1;
newpoint.nvec.p = p;
newpoint.nvec.q = q;

end
