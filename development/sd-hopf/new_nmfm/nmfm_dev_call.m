function x=nmfm_dev_call(fun,theta,deriv)
%% Evaluate history function of form sum vk t^(ak) exp(lambda_k theta)
%
% x=sum fun.v(:,i) * theta.^fun.t(i) * exp(fun.lambda(i)*theta)
%
% if present, 3rd arg deriv (0 default): differentiate so many times before
% evaluating
%
% $Id$
%%
if nargin>=3
    for i=1:deriv
        fun=nmfm_dev_deriv(fun);
    end
end
[n,nv]=size(fun.v);
nt=length(theta);
ont=ones(nt,1);
on=ones(n,1);
theta=reshape(theta,1,1,nt);
theta=theta(on,ones(nv,1),:);
v=fun.v(:,:,ont);
lam=fun.lambda(on,:,ont);
elamt=ones(n,nv,nt);
elamt(lam~=0)=exp(lam(lam~=0).*theta(lam~=0));
tp=fun.t(on,:,ont);
tp_gt_0=tp>0;
ttp=ones(n,nv,nt);
ttp(tp_gt_0)=theta(tp_gt_0).^tp(tp_gt_0);
s=v.*elamt.*ttp;
x=reshape(sum(s,2),[n,nt]);
end


