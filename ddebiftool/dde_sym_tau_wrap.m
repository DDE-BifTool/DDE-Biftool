function y=dde_sym_tau_wrap(fun,itau,order,xx,par,xdev,pdev)
%% Wrapper around automatically generated functions from symbolic differentiation
%
% $Id: dde_sym_tau_wrap.m 296 2018-09-24 21:36:56Z jansieber $
%% determine vectorized dimensions
% determine overall number of delays, which determines number of arguments
% for fun
xsize=size(xx);
nx=xsize(1);
ndelays=fun('ntau')+1;
npar=size(par,2);
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xapp=NaN(nx*(ndelays-itau),nvec);
xx=cat(1,reshape(xx(:,1:itau,:),[nx*itau,nvec]),xapp);
directional_derivative=fun('directional_derivative');
if directional_derivative
    mfrep=1;
else
    maxorder=fun('maxorder');
    mfrep=numel(xdev)/numel(xx);
end
xdev=cat(1,reshape(xdev(:,1:itau,:),[nx*itau*mfrep,nvec]),xapp);
pdev=reshape(pdev,1,npar*mfrep,[]);
if ~directional_derivative && mfrep==1
    xdev=repmat(xdev,maxorder,1);
    pdev=repmat(pdev,[1,maxorder,1]);
    mfrep=maxorder;
end
%% convert arrays to cells/argument lists in first 2 or 3 dimensions
for i=nx*ndelays:-1:1
    xxc{i}=xx(i,:);
end
for i=nx*ndelays*mfrep:-1:1
    dxxc{i}=xdev(i,:);
end
for i=npar:-1:1
    parc{i}=reshape(par(1,i,:),1,[]);
end
for i=npar*mfrep:-1:1
    dparc{i}=reshape(pdev(1,i,:),1,[]);
end
y0=fun('tau',itau,order,1,xxc{:},parc{:},dxxc{:},dparc{:});
y=NaN(1,nvec);
y(1,:)=y0;
y=reshape(y,[1,vecdim]);
end
