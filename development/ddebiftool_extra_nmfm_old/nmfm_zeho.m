function newpoint = nmfm_zeho(funcs, point, varargin)
%% Compute normal form for zero-Hopf interaction
%
%  $Id: nmfm_zeho.m 79 2015-01-02 18:42:50Z jan.sieber $
%
%%
coord = point.x;
par = point.parameter;
kind = point.kind;
n = length(coord);
ii = sqrt(-1);
newpoint = point;

if ~strcmp(kind,'zeho')
    error('NMFM_ZEHO: did not receive a zero hopf point, but %s, as argument.',kind);
end

%% Select eigenvalue pair
omega = point.omega;
if isempty(omega) || omega == 0
    warning('NMFM_ZEHO:omega',...
        'NMFM_ZEHO: omega is empty or zero, unable to compute normal form.');
    return
end

lambda0 = 0; % We could try to take the actual eigenvalue approximately equal to 0
lambda1 = ii*omega;

%% Construct characteristic matrices
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
Delta0 = ch_matrix(funcs,xx,par,0);
Delta1 = ch_matrix(funcs,xx,par,lambda1);
Delta2 = ch_matrix(funcs,xx,par,2*lambda1);
DDelta0 = ch_matrix(funcs,xx,par,lambda0,'deri',1);
DDelta1 = ch_matrix(funcs,xx,par,lambda1,'deri',1);

%% Compute nullvectors
if n == 1
    p0 = 1;
    q0 = 1;
    p1 = 1;
    q1 = 1;
else
    if nargin == 2 % use null vectors of point
        nullpoint = point;
    elseif nargin == 3 % use supplied null vectors
        nullpoint = varargin{1};
    end
    % Construct null vectors
    pp0 = []; qq0 = [];
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
        pp1 = []; 
        qq1 = [];
    end
    [p0, q0] = nmfm_border(Delta0, pp0, qq0);
    [p1, q1] = nmfm_border(Delta1, pp1, qq1);
    p0 = p0/norm(p0);
    q0 = q0/norm(q0);
    p1 = p1/norm(p1);
    q1 = q1/norm(q1);
end

%% Normalize eigenvectors
a2inv=p0*DDelta0*q0;
if a2inv<0
    q0=-q0;
    a2inv=-a2inv;
end
alpha0 = 1/sqrt(a2inv);
p0 = alpha0*p0;
q0 = alpha0*q0;

alpha1 = 1/sqrt(p1*DDelta1*q1);
p1 = alpha1*p1;
q1 = alpha1*q1;

%% Eigenfunctions
phi0 = @(theta) exp(lambda0*theta)*q0;
phi0bar = @(theta) conj(exp(lambda0*theta)*q0);
phi1 = @(theta) exp(lambda1*theta)*q1;
phi1bar = @(theta) conj(exp(lambda1*theta)*q1);

PHI0 = nmfm_handletomatrix(phi0, -taus);
PHI0BAR = nmfm_handletomatrix(phi0bar, -taus);

PHI1 = nmfm_handletomatrix(phi1, -taus);
PHI1BAR = nmfm_handletomatrix(phi1bar, -taus);

%% Some often used derivatives
D2F00 = funcs.sys_mfderi(xx,par,PHI0,PHI0);
D2F01 = funcs.sys_mfderi(xx,par,PHI0,PHI1);
D2F11b = funcs.sys_mfderi(xx,par,PHI1,PHI1BAR);

%% Normal form coefficients
%% Quadratic
g200 = (1/2)*p0*D2F00;
g110 = p1*D2F01;
g011 = p0*D2F11b;

%% Center manifold
h200 = nmfm_binv(funcs,xx,par,lambda0,q0,p0,D2F00,-p0*D2F00);
h020 = @(theta) exp(2*lambda1)*(Delta2\funcs.sys_mfderi(xx,par,PHI1,PHI1));
h110 = nmfm_binv(funcs,xx,par,lambda1,q1,p1,D2F01,-p1*D2F01);
h011 = nmfm_binv(funcs,xx,par,lambda0,q0,p0,D2F11b,-p0*D2F11b);

H200 = nmfm_handletomatrix(h200,-taus);
H020 = nmfm_handletomatrix(h020,-taus);
H110 = nmfm_handletomatrix(h110,-taus);
H011 = nmfm_handletomatrix(h011,-taus);

%% Cubic
g300 = (1/6)*p0*(3*funcs.sys_mfderi(xx,par,PHI0,H200) + funcs.sys_mfderi(xx,par,PHI0,PHI0,PHI0));
g111 = p0*(funcs.sys_mfderi(xx,par,PHI0,H011) + funcs.sys_mfderi(xx,par,PHI1BAR,H110) + ...
    funcs.sys_mfderi(xx,par,PHI1,conj(H110)) + funcs.sys_mfderi(xx,par,PHI0,PHI1,PHI1BAR));
g210 = (1/2)*p1*(funcs.sys_mfderi(xx,par,PHI1,H200) + 2*funcs.sys_mfderi(xx,par,PHI0,H110) + ...
    funcs.sys_mfderi(xx,par,PHI0,PHI0,PHI1));
g021 = (1/2)*p1*(funcs.sys_mfderi(xx,par,PHI1BAR,H020) + 2*funcs.sys_mfderi(xx,par,PHI1,H011) + ...
    funcs.sys_mfderi(xx,par,PHI1,PHI1,PHI1BAR));

%% Gavrilov
b = g200;
c = g011;
d = g110 - lambda1*g300/g200;
e = real(g210 + g110*(real(g021)/g011 - 3*g300/(2*g200) + g111/(2*g011)) - g021*g200/g011);

s = b*c;
theta0 = real(g110)/g200;

fprintf('s = %.10f, theta = %.10f.\n', s, theta0);

%% Pass on normal form coefficients
newpoint.nmfm.g200 = g200;
newpoint.nmfm.g110 = g110;
newpoint.nmfm.g011 = g011;
newpoint.nmfm.g300 = g300;
newpoint.nmfm.g111 = g111;
newpoint.nmfm.g210 = g210;
newpoint.nmfm.g021 = g021;
newpoint.nmfm.b = b;
newpoint.nmfm.c = c;
newpoint.nmfm.d = d;
newpoint.nmfm.e = e;
newpoint.nmfm.s = s;
newpoint.nmfm.theta = theta0;

newpoint.nvec.p = p1;
newpoint.nvec.q = q1;

end
