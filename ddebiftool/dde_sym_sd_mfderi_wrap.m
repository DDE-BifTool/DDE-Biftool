function y=dde_sym_sd_mfderi_wrap(fun,order,xdevlength,xx,par,xdev)
%% Wrapper around automatically generated functions from symbolic differentiation
% for higher-order directional derivatives with state-dependent delay,
% using combined functional derivatives
% $Id: dde_sym_sd_mfderi_wrap.m 296 2018-09-24 21:36:56Z jansieber $
%% determine number of arguments for fun
xx=xx(:,1);
nx=length(xx);
nf=size(xx,1);
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
out=cell(1,nf);
[out{:}]=fun('combined',1,order,nx,xxc{:},parc{:},xdevc{:});
%% The ith row gets either filled in or expanded
y=NaN(nf,1);
for i=nf:-1:1
    y(i)=out{i};
end
end
