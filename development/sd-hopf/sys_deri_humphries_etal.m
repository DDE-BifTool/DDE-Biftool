function J=sys_deri_humphries_etal(f,x,p,nx,np,v)
%% Jacobians of rhs for Humphries et al
%
% $Id: sys_deri_humphries_etal.m 20 2014-04-11 19:27:33Z jan.sieber $
%
[dim,nd,nvec]=size(x); %#ok<ASGLU>
if length(nx)==1 && isempty(np)
    if nx==0
        J=repmat(-p(5),[dim,dim,nvec]);
    elseif nx==1
        J=repmat(-p(1),[dim,dim,nvec]);
    elseif nx==2
        J=repmat(-p(2),[dim,dim,nvec]);
    end
elseif  length(np)==1 && isempty(nx)
    if np==1
        J=-x(1,2,:);
    elseif np==2
        J=-x(1,3,:);
    elseif np==5
        J=-x(1,1,:);
    else
        J=zeros([dim,1,nvec]);
    end
else
    J=df_deriv(struct('sys_rhs',f),x,p,nx,np,v);
end
end
