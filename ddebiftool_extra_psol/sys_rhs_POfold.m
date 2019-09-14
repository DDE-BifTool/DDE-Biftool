function y=sys_rhs_POfold(x,p,ip,userfuncs,sys_deri)
%% rhs of extended DDE for fold of periodic orbits
% tau contains delays ([0,p(sys_tau())]), 
% xtau_ind(i,:) contains which columns of x(1:dim,:) correspond to
% x(tau(i)+tau(1:end))
% argument sys_deri can be either numeric (say, 1e-3, then a
% finite-difference approximation is used for the linearized system), or a
% function providing the partial derivatives of the DDE's sys_rhs
%
% $Id: sys_rhs_POfold.m 369 2019-08-27 00:07:02Z jansieber $
%
dim=ip.dim;
ntau1=length(ip.orig_tau)+1;
xall=x(1:dim,:,:);
x0=xall(:,1:ntau1,:);
xd=xall(1:dim,ntau1+(1:ntau1-1),:);
v=x(dim+1:end,1:ntau1,:);
beta=p(ip.beta);
period=p(ip.period);
user_p=p(1:ip.nuserpar);
y0=userfuncs.sys_rhs(x0,user_p);
if isnumeric(sys_deri) 
    %% no user-provided derivative (sys_deri is size of deviation)
    df=@(x0,dev,ind)app_dir_deriv(@(x)userfuncs.sys_rhs(x,user_p),x0,dev,ind,sys_deri);
else
    %% sys_deri is user-provided function
    df=@(x0,dev,ind)VAopX(userfuncs.sys_deri(x0,user_p,ind-1,[],[]),dev,'*');
end
%% partial derivative wrt non-delayed term
y1=beta/period*y0+reshape(df(x0,v(:,1,:),1),size(y0));
%% add partial derivatives of all delayed terms
tau=p(ip.ext_tau);
for i=2:size(x0,2)
    deviation=v(:,i,:)+xd(:,i-1,:)*tau(i-1)*beta/period;
    y1=y1+reshape(df(x0,deviation,i),size(y0));
end
%% add partial derivatives for all parameters in nullvector
if isfield(ip,'nullparind')
    for i=1:size(ip.nullparind,1)
        i1=ip.nullparind(i,1);
        devp=p(ip.nullparind(i,2));
        y1=y1+reshape(sys_deri(x0,user_p,[],i1,[])*devp,size(y1));
    end
end
%% assemble overall r.h.s.
y=cat(1,y0,y1);
end
