function y=dde_sym_wrap(fun,action,ind,order,nf,xx,par,xdev,pdev)
%% Wrapper around automatically generated functions from symbolic differentiation
%
% $Id: dde_sym_wrap.m 172 2017-03-06 10:26:06Z jansieber $
%% determine vectorized dimensions
xsize=size(xx);
ndelays=xsize(2);
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xx=reshape(xx,[xsize(1),ndelays,nvec]);
%% convert arrays to cells/argument lists in first 2 dimensions
xxc=num2cell(xx,3);
parc=num2cell(par,3);
dxxc=num2cell(xdev,3);
dparc=num2cell(pdev,3);
out=cell(1,nf);
[out{:}]=fun(action,ind,order,nf,xxc{:},parc{:},dxxc{:},dparc{:});
%% The ith row gets either filled in or expanded
y=NaN(nf,nvec);
for i=nf:-1:1
    y(i,:)=out{i};
end
y=reshape(y,[nf,vecdim]);
end
