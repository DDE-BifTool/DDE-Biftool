function J=sys_dtau_humphries_etal(sys_tau,itau,x,p,nx,np)
%% Jacobians of delays for Humphries et al
%
% $Id: sys_dtau_humphries_etal.m 20 2014-04-11 19:27:33Z jan.sieber $
%
[dim,nd,nvec]=size(x); %#ok<ASGLU>
if length(nx)==1 && isempty(np)
    if nx==0
        J=repmat(p(6),[1,1,nvec]);
    else
        J=zeros([1,1,nvec]);
    end
elseif length(np)==1 && isempty(nx)
    if np==6
        J=x(1,1,:);
    elseif np-itau==2
        J=ones([1,1,nvec]);
    else
        J=zeros([1,1,nvec]);
    end
else
    J=df_derit(struct('sys_tau',sys_tau),itau,x,p,nx,np);
end
end
