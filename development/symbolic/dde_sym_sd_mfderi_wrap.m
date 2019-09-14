function y=dde_sym_sd_mfderi_wrap(fun,order,xdevlength,xx,par,xdev)
%% Wrapper around automatically generated functions from symbolic differentiation
% for higher-order directional derivatives with state-dependent delay,
% using combined functional derivatives
% $Id$
%% determine number of arguments for fun
xx=xx(:,1);
nx=length(xx);
xxc=num2cell(xx,2);
parc=num2cell(par,3);
xdevc=cell(1,xdevlength);
getk=@(x,k)x(k);
for k=1:nx
    xxc{k}=@(t)xx(k);
    for i=1:order
        xdevc{(i-1)*nx+k}=@(t)getk(xdev(i-1,t),k);
    end
end
y=fun('combined',1,order,nx,xxc{:},parc{:},xdevc{:});
end
