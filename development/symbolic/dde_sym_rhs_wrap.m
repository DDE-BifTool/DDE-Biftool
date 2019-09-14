function y=dde_sym_rhs_wrap(fun,order,xx,par,xdev,pdev)
%% Wrapper around automatically generated functions from symbolic differentiation
%
% $Id: dde_sym_wrap.m 168 2017-03-03 22:04:32Z jansieber $
%% determine vectorized dimensions
xsize=size(xx);
nf=xsize(1);
ndelays=xsize(2);
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xx=reshape(xx,[nf*ndelays,nvec]);
xdev=reshape(xdev,[nf*ndelays,nvec]);
%% convert arrays to cells/argument lists in first 2 dimensions
for i=nf*ndelays:-1:1
    xxc{i}=xx(i,:);
    dxxc{i}=xdev(i,:);
end
for i=length(par):-1:1
    parc{i}=par(i);
    dparc{i}=pdev(i);
end
out=cell(1,nf);
[out{:}]=fun('rhs',1,order,nf,xxc{:},parc{:},dxxc{:},dparc{:});
%% The ith row gets either filled in or expanded
y=NaN(nf,nvec);
for i=nf:-1:1
    y(i,:)=out{i};
end
y=reshape(y,[nf,vecdim]);
end
