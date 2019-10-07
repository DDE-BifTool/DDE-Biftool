function [r,J]=dde_hcli_adjoint_proj(funcs,hcli,free_par,varargin)
%% computes projection of end of collocation segment with adjoint
% at fixed point hcli.x2
%
% Note that the vectors wj are the complex conjugates of the
% adjoint eigenvectors such that w.'*Delta(lambda)=0
%%
default={'time',1};
options=dde_set_options(default,varargin,'pass_on');
[ind,len]=dde_ind_from_point(hcli,free_par); %#ok<ASGLU>
n2=length(hcli.lambda_w);
x2=hcli.x2;
dim=length(x2);
lambda_w=hcli.lambda_w;
vec_w=hcli.w;
xt=hcli.profile;
stst=dde_stst_create('x',x2,'parameter',hcli.parameter);
Fvals=dde_stst_sysvals(funcs,stst,free_par,...
    'output',{'dfdx','dfdp','dfdxx','dfdxp','dtaudx','dtaudp'});
tau=Fvals.tau(2:end);
ntau=length(tau);
dfdx=Fvals.dfdx(:,:,2:end);
dfdx=reshape(dfdx,dim,dim*ntau);
dfdx_xp=cat(4,Fvals.dfdxx(:,:,2:end,:),Fvals.dfdxp(:,:,2:end,:));
dfdx_xp=permute(reshape(dfdx_xp,dim,dim*ntau,[]),[1,3,2]);
wdfdx_xp=reshape(vec_w.'*reshape(dfdx_xp,dim,[]),[],dim*ntau);
dtau_xp=cat(2,Fvals.dtaudx(2:end,:),Fvals.dtaudp(2:end,:));
lhs=funcs.lhs_matrix(dim);
Jp=repmat(p_axpy(0,hcli,[]),n2,1);
for i=n2:-1:1
    lam=lambda_w(i);
    w=vec_w(:,i);
    [int,J_int,ip]=dde_coll_intexp(hcli,tau,lam,options.time);
    Jcw=lhs*(xt(:,end)-x2)+dfdx*int;
    % residual
    cr(i)=w.'*Jcw;
    r(:,i)=[real(cr(i));imag(cr(i))];
    % complex jacobian
    Jp(i).w(:,i)=Jcw;
    Jp(i).profile=reshape(w.'*dfdx*J_int(:,ip.profile),size(hcli.profile));
    Jp(i).profile(:,end)=Jp(i).profile(:,end)+lhs.'*w;
    Jp(i).x2=-lhs.'*w;
    Jp(i).lambda_w(i)=w.'*dfdx*J_int(:,ip.lambda);
    Jp(i).period=w.'*dfdx*J_int(:,ip.period);
    int_dxp=wdfdx_xp(i,:)*int+...
        w.'*dfdx*J_int(:,ip.tau)*dtau_xp;
    Jp(i).x2=Jp(i).x2+int_dxp(1:dim);
    Jp(i).parameter(free_par)=int_dxp(dim+(1:length(free_par)));
end
r=r(:);
Jc=arrayfun(@(p)dde_complex_from_point(p,free_par),Jp,'uniformoutput',false);
J=cat(2,Jc{:}).';
end