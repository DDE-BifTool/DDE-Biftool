function y=dde_sym_rhs_wrap(fun,order,xx,par,xdev,pdev)
%% Wrapper around automatically generated functions from symbolic differentiation
%
% $Id: dde_sym_rhs_wrap.m 315 2019-01-29 19:42:21Z jansieber $
%% determine vectorized dimensions
xsize=size(xx);
nf=xsize(1);
ndelays=xsize(2);
npar=size(par,2);
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xx=reshape(xx,[nf*ndelays,nvec]);
mfrep=numel(xdev)/numel(xx);
xdev=reshape(xdev,[nf*ndelays*mfrep,nvec]);
pdev=reshape(pdev,1,npar*mfrep,[]);
%% convert arrays to cells/argument lists in first 2 or 3 dimensions
for i=nf*ndelays:-1:1
    xxc{i}=xx(i,:);
end
for i=nf*ndelays*mfrep:-1:1
    dxxc{i}=xdev(i,:);
end
for i=npar:-1:1
    parc{i}=reshape(par(1,i,:),1,[]);
end
for i=npar*mfrep:-1:1
    dparc{i}=reshape(pdev(1,i,:),1,[]);
end
out=cell(1,nf);
if ~isempty(xdev)
    %% if multi-derivative is called as directional derivative, replicate deviations
    if order>mfrep && ~fun('directional_derivative')
        dxxc=repmat(dxxc,1,order/mfrep);
        dparc=repmat(dparc,1,order/mfrep);
    end
    %% call derivative
    [out{:}]=fun('rhs',1,order,nf,xxc{:},parc{:},dxxc{:},dparc{:});
else
    %% call function
    [out{:}]=fun('rhs',1,order,nf,xxc{:},parc{:});
end
%% The ith row gets either filled in or expanded
y=NaN(nf,nvec);
for i=nf:-1:1
    y(i,:)=out{i};
end
y=reshape(y,[nf,vecdim]);
end
