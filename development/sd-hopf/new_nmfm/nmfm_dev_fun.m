function fun=nmfm_dev_fun(v,varargin)
%% Create history function of form sum vk t^(ak) exp(lambda_k theta)
%
% x=sum fun.v(:,i) * theta.^fun.t(i) * exp(fun.lambda(i)*theta)
%
% $Id$
%% Optional arguments
default={'lambda',0,'t',0};
options=dde_set_options(default,varargin,'pass_on');
nv=size(v,2);
lambda=options.lambda(:).';
if length(lambda)~=nv
    error('nmfm_dev_fun:args',...
        'nmfm_dev_fun: length of lambda must equal number of v''s');
end
t=options.t(:)';
if isempty(t)
    t=zeros(size(lambda));
elseif length(t)~=nv
    error('nmfm_dev_fun:args',...
        'nmfm_dev_fun: length of t must equal number of v''s');
end
%% combine duplicates (if lambda_k and t_k are exactly equal)
[tl,~,ic]=unique([t;lambda].','rows');
tnew=tl(:,1).';
lambdanew=tl(:,2).';
vnew=zeros(size(v,1),length(lambdanew));
for i=1:length(ic)
    vnew(:,ic(i))=vnew(:,ic(i))+v(:,i);
end
fun=struct('v',vnew,'lambda',lambdanew,'t',tnew);
end

    
    