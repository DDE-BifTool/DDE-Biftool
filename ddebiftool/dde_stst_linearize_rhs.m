function [A,taus]=dde_stst_linearize_rhs(funcs,point,varargin)
%% linearize r.h.s. in stst
%
% $Id: dde_stst_linearize_rhs.m 302 2018-09-28 18:06:58Z jansieber $
%%
default={'free_par','none'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
taus=dde_stst_delays(funcs,point,pass_on{:});
m=length(taus);
xx=point.x(:,ones(m,1));
d_x=ischar(options.free_par);
if d_x % 1st or 2nd  deriv wrt x
    for k=m:-1:1
        A(:,:,k)=funcs.sys_deri(xx,point.parameter,k-1,[],[]);
    end
else % d/dp or d^2/(dx dp)
    nf=length(options.free_par);
    A=[];
    for k=nf:-1:1
        A(:,k)=funcs.sys_deri(xx,point.parameter,[],options.free_par(k),[]);
    end
end
end
