function y=dde_sym_tau_wrap(fun,itau,order,xx,par,xdev,pdev)
%% Wrapper around automatically generated functions from symbolic differentiation
%
% $Id: dde_sym_wrap.m 168 2017-03-03 22:04:32Z jansieber $
%% determine vectorized dimensions
% determine overall number of delays, which determines number of arguments
% for fun
xsize=size(xx);
nx=xsize(1);
ndelays=fun('ntau')+1;
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xapp=NaN(nx*(ndelays-itau),nvec);
xx=cat(1,reshape(xx(:,1:itau,:),[nx*itau,nvec]),xapp);
xdev=cat(1,reshape(xdev(:,1:itau,:),[nx*itau,nvec]),xapp);
%% convert arrays to cells/argument lists in first 2 dimensions
for i=nx*ndelays:-1:1
    xxc{i}=xx(i,:);
    dxxc{i}=xdev(i,:);
end
for i=length(par):-1:1
    parc{i}=par(i);
    dparc{i}=pdev(i);
end
y=fun('tau',itau,order,1,xxc{:},parc{:},dxxc{:},dparc{:});
y=reshape(y,[1,vecdim]);
end
