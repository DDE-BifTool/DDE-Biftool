function [J,Jd,ind]=dde_coll_intexp(coll,tau,lambda,nc)
%% compute int_{-tau/T}^0 T exp(-lambda(T*s+tau))*u(1+s) ds
% needed for spectral projection of s->u(1+s) and its derivatives
%
% inputs
% * coll: collocation structure for u (with period)
% * tau: delay
% * lambda: exponent
% * nc: number of Jacobian columns
% 
%
% output:
%
% * J: integral
% * Jd: derivative, (rows are rows of coll.profile, 
% * ind: indices for columns in derivatives are
% ind.lambda,ind.profile,ind.tau,ind.period,
%
% $Id: dde_coll_intexp.m 375 2019-09-14 14:26:50Z jansieber $
%%
[dim,nt]=size(coll.profile);
ntau=length(tau);
ind=struct('lambda',1,'period',2,'profile',reshape(2+(1:dim*nt),dim,nt),...
    'tau',2+dim*nt+(1:ntau));
Jc=cell(ntau,1);
Jdc=cell(ntau,1);
for i=1:ntau
    [Jc{i},Jdc{i}]=intexp_loc(coll,tau(i),lambda,...
        setfield(ind,'tau',ind.tau(i)),nc); %#ok<*SFLD>
end
J=cat(1,Jc{:});
Jd=cat(1,Jdc{:});
end
%%
function [J,Jd]=intexp_loc(coll,tau,lambda,ind,nc)
assert(tau<coll.period);
[dim,nt]=size(coll.profile);
s=coll.mesh-1;
T=coll.period;
explam=exp(-lambda*(T*s+tau));
%% explam could be very large for initial times, while tau could be small
% so, avoid 0*inf by setting explam to zero for all values
explam(T*s+tau<0)=0;
c0=setfield(coll,'profile',zeros(1,nt));
utau=dde_coll_eva(coll.profile,coll.mesh,1-tau/T,coll.degree);
[dum,W]=dde_coll_profile_dot(c0,c0,'bd',[1-tau/T,1]); %#ok<ASGLU>
J=T*coll.profile*W'*explam(:);
Jprof=kron(T*explam(:).'*W,speye(dim));
sc=@(y)coll.profile*W'*y(:);
JT=sc((1-lambda*T*s).*explam)-tau/T*utau;
Jtau=sc(-lambda*T*explam)+utau;
Jlambda=sc(-T*(T*s+tau).*explam);
ic=repmat([ind.lambda,ind.tau,ind.period,ind.profile(:)'],dim,1);
ir=repmat((1:dim)',1,size(ic,2));
vals=     [   Jlambda,   Jtau,     JT,       Jprof];
Jd=sparse(ir(:),ic(:),vals(:),dim,nc);
end
