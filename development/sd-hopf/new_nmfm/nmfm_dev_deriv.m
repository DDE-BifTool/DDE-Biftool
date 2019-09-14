function dfun=nmfm_dev_deriv(fun)
%% Differentiate history function (once)
%
% $Id$
%%
n=size(fun.v,1);
lambda=fun.lambda;
v=fun.v.*lambda(ones(n,1),:);
t=fun.t;
t_gt_0=t>0;
if ~any(t_gt_0)
    dfun=nmfm_dev_fun(v,'lambda',lambda,'t',t);
    return
end
tex=t(t_gt_0);
vex=fun.v(:,t_gt_0).*tex(ones(n,1),:);
lex=lambda(:,t_gt_0);
tex=tex-1;
v=[v,vex];
lambda=[lambda,lex];
t=[t,tex];
dfun=nmfm_dev_fun(v,'lambda',lambda,'t',t);
end

