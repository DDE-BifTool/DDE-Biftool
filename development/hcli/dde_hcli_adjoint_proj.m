function [r,J]=dde_hcli_adjoint_proj(Fvals,hcli,free_par)
[ind,len]=dde_ind_from_point(hcli,free_par);
n2=length(hcli.lambda_w);
x2=hcli.x2;
dim=length(x2);
lambda_w=hcli.lambda_w;
vec_w=hcli.w;
xt=hcli.profile;
T=hcli.period;
tau=Fvals.tau2(2:end);
ntau=length(tau);
dfdx=Fvals.dfdx2(:,:,2:end);
dfdx=reshape(dfdx,dim,dim*ntau);
dfdx_xp=cat(4,Fvals.dfdxx2(:,:,2:end,:),Fvals.dfdxp2(:,:,2:end,:));
dfdx_xp=permute(reshape(dfdx_xp,dim,dim*ntau,[]),[1,3,2]);
dtau_xp=cat(2,Fvals.dtaudx2(2:end,:),Fvals.dtaudp2(2:end,:));
lhs=Fvals.lhs_mat;
Jp=repmat(p_axpy(0,hcli,[]),n2,1);
for i=n2:-1:1
    lam=lambda_w(i);
    w=vec_w(:,i);
    [int,J_int,ip]=dde_coll_intexp(hcli,tau,lam,len);
    Jcw=lhs*(xt(:,end)-x2)+dfdx*int;
    % residual
    cr(i)=w.'*Jcw;
    r(:,i)=[real(cr(i));imag(cr(i))];
    % complex jacobian
    Jp(i).w(:,i)=Jcw;
    Jp(i).profile=w.'*dfdx*J_int(:,ip.profile);
    Jp(i).profile(:,end)=Jp(i).profile(:,end)+lhs.'*w;
    Jp(i).x2=-lhs.'*w;
    Jp(i).lambda_w(i)=w.'*dfdx*J_int(:,ip.lambda);
    Jp(i).period=w.'*dfdx*J_int(:,ip.period);
    Jp(i).x2=Jp(i).x2+...
        reshape(w.'*reshape(dfdx_xp,dim,[]),dim,dim*ntau)*int+...
        w.'*dfdx*J_int(:,ip.tau)*dtau_xp;
end
end