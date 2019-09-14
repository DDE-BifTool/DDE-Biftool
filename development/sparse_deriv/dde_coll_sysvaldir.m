function [vals,nf]=dde_coll_sysvaldir(order,funcs,xx,p,yy,q,shape)
%% evaluate directional derivatives of rhs in a vectorized form
% xx has shape dim nc d x 1, yy has shape dim x nc x d x ndev, p has
% shape 1 x np, q has shape 1 x np x ndev
%
% output is shape dim x nc x ndev
xx=reshape(xx,shape);
dim=shape(1);
nc=shape(2);
d=shape(3);
yy=reshape(full(yy),dim,nc,d,[]);
ndev=size(yy,4);
xx=permute(xx,[1,3,2]);
if order==0
    vals=funcs.sys_rhs(xx,p);
    nf=size(vals,1);
    vals=reshape(vals,nf*nc,1);
    return
end
q=repmat(reshape(full(q),1,[],1,ndev),1,1,nc,1);
xx=repmat(xx,1,1,1,ndev);
p=repmat(p,1,1,nc,ndev);
p=reshape(p,1,[],nc*ndev);
q=reshape(q,1,[],nc*ndev);
xx=reshape(xx,dim,d,nc*ndev);
yy=permute(yy,[1,3,2,4]);
yy=reshape(yy,dim,d,nc*ndev);
vals=funcs.sys_dirderi{order}(xx,p,yy,q);
nf=size(vals,1);
vals=reshape(vals,[],ndev);
end