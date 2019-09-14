function bt = nmfm_bt(funcs, bt, varargin)

if ~strcmp(bt.kind,'BT')
   display(bt.kind);
   error('NMFM_BT: did not receive a bt point as argument.');
end

sys_tau=funcs.sys_tau;
sys_mfderi=funcs.sys_mfderi;

x=bt.x;
par=bt.parameter;
tau=[0 par(sys_tau())];
m=length(tau)-1;
xx=x(:,ones(m+1,1));

% n=length(x);
% 
% l=0;
% D=l*eye(n);
% dD=eye(n);
% ddD=zeros(n);
% dddD=zeros(n);
% 
% for j=0:m
%   B=sys_deri(xx,par,j,[],[])*exp(-l*tau(j+1));
%   D=D-B;
%   dD=dD+tau(j+1)*B;
%   ddD=ddD-tau(j+1)^2*B;
%   dddD=dddD+tau(j+1)^3*B;
% end

D=ch_matrix(funcs,bt.x,bt.parameter,0);
dD=ch_matrix(funcs,bt.x,bt.parameter,0,'deri',1);
ddD=ch_matrix(funcs,bt.x,bt.parameter,0,'deri',2);
dddD=ch_matrix(funcs,bt.x,bt.parameter,0,'deri',3);

%% normalize the jordan chains vectors
q0=bt.q0; q1=bt.q1; p0=bt.p0; p1=bt.p1;

% check conditions
% D*bt.q0
% dD*bt.q0+D*bt.q1
% p1'*D
% p1'*dD+p0'*D

alpha=p0'*dD*q0+1/2*p1'*ddD*q0;
p0=p0/alpha;
p1=p1/alpha;

% beta=p1'*dD*q1+1/2*p1'*ddD*q0;
% p1=p1/beta;

delta=-(p0'*dD*q1+1/2*p0'*ddD*q0+1/2*p1'*ddD*q1+1/6*p1'*dddD*q0);
q1=q1+delta*q0;

% check
% p0'*dD*q0+1/2*p1'*ddD*q0
% p1'*dD*q1+1/2*p1'*ddD*q0
% bt.p1'*dD*bt.q0
% p0'*dD*q1+1/2*p0'*ddD*q0+1/2*p1'*ddD*q1+1/6*p1'*dddD*q0

%% calculate normal form coefficients
phi_0 = @(theta) q0;
phi_1 = @(theta) theta*q0+q1;

PHI_0 = nmfm_handletomatrix(phi_0, -tau);
PHI_1 = nmfm_handletomatrix(phi_1, -tau);

bt.nmfm.a2=1/2*p1'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0);
bt.nmfm.b2=p0'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0)+...
    p1'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_1);

% bt.nmfm.a3=1/6*p1'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_0);
% bt.nmfm.b3=1/2*p1'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_1) ...
%     + 1/2*p0'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_0);


end