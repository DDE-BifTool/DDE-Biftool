function y=sys_rhs_MWFold(x,p,ind,userfuncs,sys_deri)
%% rhs of extended DDE for fold of modulated waves
%
% tau contains delays ([0,p(sys_tau())]), 
% xtau_ind(i,:) contains which columns of x(1:dim,:) correspond to
% x(tau(i)+tau(1:end))
% beta is derivative wrt period
% rho is rotation frequency
%
% $Id: sys_rhs_MWFold.m 369 2019-08-27 00:07:02Z jansieber $
%
y=sys_rhs_POfold(x,p,ip,userfuncs,sys_deri);
%% partial derivative wrt rho
x0=x(1:ind.dim,:);
devx=0*x0;
devp=zeros(1,ind.nuserpar);
for i=1:size(ip.nullparnames,1)
    i1=ind.(ip.nullparnames{i,1});
    i2=ind.(ip.nullparnames{i,2});
    devp(1,i1)=p.parameter(i2);
end
df=userfuncs.sys_dirderi{1}(x0,p,devx,devp);
y(ind.dim+(1:ind:dim))=y(ind.dim+(1:ind:dim))+df;
end
