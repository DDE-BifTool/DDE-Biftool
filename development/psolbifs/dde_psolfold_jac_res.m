function [J,res]=dde_psolfold_jac_res(funcs,pt,free_par,method,varargin)
default={'hdev',1e-4};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
extfuncs=funcs;
if ~funcs.tp_del
    itau=length(funcs.sys_tau());
    extfuncs.sys_tau=@()[itau,itau];
    npar=length(pt.parameter);
    npar_var=length(pt.parameter_var);
    dtaudp=
    extfuncs.sys_rhs=@(x,p)sys_rhs_psolfold(x,p(1:npar),p(npar+(1:npar_var)),...
        funcs,free_par(1:end-1),dtaudp,options.hdev);
        
else
end
end
function y=sys_rhs_psolfold(x,p,phat,userfuncs,free_par,dtaudp,numderi)
%% rhs of extended DDE for fold of periodic orbits
% first x(:,2:ntau+1,:) contain x(t-tau(j)), ntau+(1:ntau)
% contain x'(t-tau(j)), j=1:ntau
% argument sys_deri can be either numeric (say, 1e-3, then a
% finite-difference approximation is used for the linearized system), or a
% function providing the partial derivatives of the DDE's sys_rhs
%
% $Id: dde_psolfold_jac_res.m 340 2019-05-09 19:50:29Z jansieber $
%
itau=userfuncs.sys_tau();
dim=size(x,1)/2;
ntau=length(itau)+1;
tau=[0,p(itau)]';
xall=x(1:dim,:,:);
x0=xall(:,1:ntau,:);            % x(t-tau)
v=x(dim+1:end,1:ntau,:);        % v(t-tau)
xp=xall(1:dim,ntau+(1:ntau),:); % x'(t-tau)
y0=funcs.sys_rhs(x0,p);
ry=@(y)reshape(y,size(y0));
beta=phat(1);
devpars=reshape(phat(2:end),[],1);
devtau=dtaudp*devpars;
if isnumeric(numderi) 
    %% no user-provided derivative (numderi is size of deviation)
    dfx=@(dev,ind)ry(app_dir_deriv(@(xa)funcs.sys_rhs(xa,p),x0,dev,ind,numderi));
    dfp=@(dev,ind)ry(app_dir_deriv(@(pa)funcs.sys_rhs(x,pa),p ,dev,ind,numderi));
else
    %% sys_deri is user-provided function
    dfx=@(dev,ind)ry(VAopX(funcs.sys_deri(x,p,ind-1,[],[]),dev,'*'));
    dfp=@(dev,ind)ry(funcs.sys_deri(x,p,[],ind,[])*dev);
end
%% partial derivative wrt period
y1=beta*xp(:,1,:);
%% add partial derivatives of all delayed terms
for i=1:size(x0,2)
    fdev=xp(:,i,:)*(tau(i)*beta-devtau(i));
    deviation=v(:,i,:)+fdev;
    y1=y1+dfx(deviation,i);
end
for i=1:length(devpars)
    y1=y1+dfp(devpars(i),free_par(i));
end
y=cat(1,y0,y1);
end
