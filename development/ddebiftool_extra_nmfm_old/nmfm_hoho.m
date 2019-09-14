function newpoint = nmfm_hoho(funcs, point, varargin)
%% Compute normal form of Double Hopf point
%
% $Id: nmfm_hoho.m 151 2017-02-14 11:38:46Z jansieber $
%
%%
coord = point.x;
par = point.parameter;
kind = point.kind;
n = length(coord);
ii = sqrt(-1);
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

lambda1 = ii*omega1;
lambda2 = ii*omega2;

%% Get delays
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

% Compute nullvectors
if n == 1
   p1 = 1;
   q1 = 1;
   p2 = 1;
   q2 = 1;
else
   if nargin == 2 % use null vectors of point
      nullpoint = point;
   elseif nargin == 3 % use supplied null vectors
      nullpoint = varargin{1};
   end
   % Construct null vectors
   % Construct null vectors
   pp2 = []; qq2 = [];
   if isfield(nullpoint, 'nvec')
      if isfield(nullpoint.nvec,'p') 
          pp1 = nullpoint.nvec.p; 
      else
          pp1 = [];
      end
      if isfield(nullpoint.nvec,'q')
          qq1 = nullpoint.nvec.q; 
      else
          qq1 = []; 
      end
   else
      pp1 = []; qq1 = [];
   end
   [p1, q1] = nmfm_border(ch_matrix(funcs,xx,par,lambda1), pp1, qq1);
   [p2, q2] = nmfm_border(ch_matrix(funcs,xx,par,lambda2), pp2, qq2);
   p1 = p1/norm(p1);
   q1 = q1/norm(q1);
   p2 = p2/norm(p2);
   q2 = q2/norm(q2);
end

if isempty(p2) || isempty(q2) || isempty(q1) || isempty(p1)
   fprintf('NMFM_HOHO: could not construct null vectors.\n');
   return;
end

% Normalize eigenvectors
alpha1 = 1/sqrt(p1*ch_matrix(funcs,xx,par,lambda1,'deri',1)*q1);
p1 = alpha1*p1;
q1 = alpha1*q1;

alpha2 = 1/sqrt(p2*ch_matrix(funcs,xx,par,lambda2,'deri',1)*q2);
p2 = alpha2*p2;
q2 = alpha2*q2;

% Eigenfunctions
phi1 = @(theta) exp(lambda1*theta)*q1;
phi1bar = @(theta) conj(exp(lambda1*theta)*q1);

phi2 = @(theta) exp(lambda2*theta)*q2;
phi2bar = @(theta) conj(exp(lambda2*theta)*q2);

PHI1 = nmfm_handletomatrix(phi1, -taus);
PHI1BAR = nmfm_handletomatrix(phi1bar, -taus);

PHI2 = nmfm_handletomatrix(phi2, -taus);
PHI2BAR = nmfm_handletomatrix(phi2bar, -taus);

% Normal form coefficients
% Quadratic center manifold
h1100 = @(theta) ch_matrix(funcs,xx,par,0)\funcs.sys_mfderi(xx,par,PHI1,PHI1BAR);
h2000 = @(theta) exp(2*lambda1*theta)*(ch_matrix(funcs,xx,par,2*lambda1)\funcs.sys_mfderi(xx,par,PHI1,PHI1));
h1010 = @(theta) exp((lambda1+lambda2)*theta)*(ch_matrix(funcs,xx,par,lambda1+lambda2)\ ...
   funcs.sys_mfderi(xx,par,PHI1,PHI2));
h1001 = @(theta) exp((lambda1-lambda2)*theta)*(ch_matrix(funcs,xx,par,lambda1-lambda2)\ ...
   funcs.sys_mfderi(xx,par,PHI1,PHI2BAR));
h0020 = @(theta) exp(2*lambda2*theta)*(ch_matrix(funcs,xx,par,2*lambda2)\funcs.sys_mfderi(xx,par,PHI2,PHI2));
h0011 = @(theta) ch_matrix(funcs,xx,par,0)\funcs.sys_mfderi(xx,par,PHI2,PHI2BAR);

H1100 = nmfm_handletomatrix(h1100,-taus);
H2000 = nmfm_handletomatrix(h2000,-taus);
H1010 = nmfm_handletomatrix(h1010,-taus);
H1001 = nmfm_handletomatrix(h1001,-taus);
H0020 = nmfm_handletomatrix(h0020,-taus);
H0011 = nmfm_handletomatrix(h0011,-taus);

% Cubic normal form
g2100 = (1/2)*p1*(2*funcs.sys_mfderi(xx,par,H1100,PHI1) + funcs.sys_mfderi(xx,par, H2000, PHI1BAR) + ...
   funcs.sys_mfderi(xx,par,PHI1,PHI1,PHI1BAR));
g1011 = p1*(funcs.sys_mfderi(xx,par,H0011,PHI1) + funcs.sys_mfderi(xx,par,H1001,PHI2) + ...
   funcs.sys_mfderi(xx,par,H1010,PHI2BAR) + funcs.sys_mfderi(xx,par, PHI1,PHI2,PHI2BAR));
g1110 = p2*(funcs.sys_mfderi(xx,par,conj(H1001),PHI1) + funcs.sys_mfderi(xx,par, H1010,PHI1BAR) + ...
   funcs.sys_mfderi(xx,par, H1100, PHI2) + funcs.sys_mfderi(xx,par, PHI1, PHI1BAR,PHI2));
g0021 = (1/2)*p2*(2*funcs.sys_mfderi(xx,par,H0011,PHI2) + funcs.sys_mfderi(xx,par, H0020,PHI2BAR) + ...
   funcs.sys_mfderi(xx,par,PHI2,PHI2,PHI2BAR));

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