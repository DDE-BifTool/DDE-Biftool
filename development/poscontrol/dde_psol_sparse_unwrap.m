%% convert structs to sparse matrices for psol
%
% $Id$
%%
function [W,mesh]=dde_psol_sparse_unwrap(W,basemesh,degree)
% W{t_i} are structs, check for minimal index and shift to create
% sparse matrices
imin=1;
tm=length(basemesh)-1;
xmin=0;
xmax=1;
for t_i=1:length(W)
    imin=min(W{t_i}.imin,imin);
    xmin=min([W{t_i}.x,xmin]);
    xmax=max([W{t_i}.x,xmax]);
end
%% extend mesh
xfmin=floor(xmin);
xfmax=floor(xmax);
mesh=repmat(basemesh(:),1,xfmax-xfmin+1)+repmat(xfmin:xfmax,tm+1,1);
c_is_int=mesh==round(mesh);
mesh(1,c_is_int(1,:)&mesh(1,:)>=1)=NaN;
mesh(end,c_is_int(end,:)&mesh(end,:)<=0)=NaN;
mesh=mesh(:)';
mesh=mesh(~isnan(mesh));
coarsemesh=mesh(1:degree:end);
icmin=find(coarsemesh<=xmin,1,'last');
icmax=find(coarsemesh>=xmax,1,'first');
mesh=mesh((icmin-1)*degree+1:(icmax-1)*degree+1);
%% shift matrix indices and create sparse matrices
for t_i=1:length(W)
    W{t_i}.ic=W{t_i}.ic-imin+1;
    W{t_i}=sparse(W{t_i}.ir,W{t_i}.ic,W{t_i}.vals,...
            W{t_i}.nr,W{t_i}.nc-imin+1);
end
end
