function bt = nmfm_bt(funcs, bt, varargin)
%% Bogdanov-Takens normal form coefficients
%
% This function assumes that the nullvectors of bt have been computed
% already.
%
% $Id: nmfm_bt.m 309 2018-10-28 19:02:42Z jansieber $
%%
default={'debug',false};
options=dde_set_options(default,varargin,'pass_on');
if ~strcmp(bt.kind,'BT')
   display(bt.kind);
   error('NMFM_BT: did not receive a bt point as argument.');
end
D=ch_matrix(funcs,bt.x,bt.parameter,0);
dD=ch_matrix(funcs,bt.x,bt.parameter,0,'deri',1);
ddD=ch_matrix(funcs,bt.x,bt.parameter,0,'deri',2);
dddD=ch_matrix(funcs,bt.x,bt.parameter,0,'deri',3);

%% normalize the Jordan chains vectors
% q0'q0=1 and q0'*q1=0 from defining system (that is, q1 could be 0)
q0=bt.nvec.q0; 
q1=bt.nvec.q1; 
p0=bt.nvec.p0; 
p1=bt.nvec.p1;

beta= -( p0'*dD*q1 + 1/2*p0'*ddD*q0 + 1/2*p1'*ddD*q1 + 1/6*p1'*dddD*q0 ) /...
    (p1'*dD*q1+1/2*p1'*ddD*q0);
p0=p0+beta*p1;
alpha =  p0'*dD*q0 +  1/2 * p1'*ddD*q0;
p0=p0/alpha;
p1=p1/alpha;
%% Check conditions and normalizations
% note that q1 or p0 may be 0.
if options.debug
    small=@(x)norm(x)<1e-12;
    assert(small(D*q0));          %  A phi0=0
    assert(small(dD*q0+D*q1));    %  A phi1+phi0=0
    assert(small(p1'*D));         % phisun1 A = 0
    assert(small(p0'*D+p1'*dD));  % phisun1 + phisun0 A = 0
    assert(small(q0'*q0-1));      % |phi0(0)|=1 (this is arbitrary)
    assert(small(q0'*q1));        % phi1(0)^T phi0(0)=0 (this is arbitrary)
    assert(small(p0'*dD*q0 + 1/2*p1'*ddD*q0 - 1)); % <phisun0,phi0> = 1
    assert(small(p1'*dD*q1 + 1/2*p1'*ddD*q0 - 1)); % <phisun1,phi1> = 1
    assert(small(p0'*dD*q1 + 1/2*p0'*ddD*q0 + 1/2*p1'*ddD*q1 + 1/6*p1'*dddD*q0)); % <phisun0,phi1> = 0
end
%% abbreviate derivative
F=nmfm_deriv_define(funcs,bt,'free_pars',[],varargin{:});
par0=bt.parameter(:)*0;
dev0=@(v)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))]);
devl=@(v,lambda)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',lambda); %#ok<NASGU>
devlt=@(v,lambda,t)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',lambda,'t',t);

% delta=-(p0'*dD*q1+1/2*p0'*ddD*q0+1/2*p1'*ddD*q1+1/6*p1'*dddD*q0);
% q1=q1+delta*q0;

% check
% p0'*dD*q0+1/2*p1'*ddD*q0
% p1'*dD*q1+1/2*p1'*ddD*q0
% bt.p1'*dD*bt.q0
% p0'*dD*q1+1/2*p0'*ddD*q0+1/2*p1'*ddD*q1+1/6*p1'*dddD*q0

%% calculate normal form coefficients
phi_0=dev0(q0);
phi_1=devlt([q1,q0],[0,0],[0,1]);
Bphi02=F.B(phi_0,phi_0);
bt.nmfm.a2=1/2*p1'*Bphi02;
bt.nmfm.b2= p0'*Bphi02 + p1'*F.B(phi_0,phi_1);
bt.nvec.p1=p1;
bt.nvec.p0=p0;
% bt.nmfm.a3=1/6*p1'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_0);
% bt.nmfm.b3=1/2*p1'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_1) ...
%     + 1/2*p0'*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_0);


end