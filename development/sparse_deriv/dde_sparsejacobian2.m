function df=dde_sparsejacobian2(f,x,varargin)
%% compute Jacobian using finite differences (optionally prescribed pattern)
% optional inputs: 
%
% * 'h': deviation (default 1e-5)
% * 'order': order of finite difference formulas (1 or 2, default 2)
% * 'f0': f(x) (if f(x) has already been computed, default [] to recompute)
% * 'vectorised', if f can be called with several vectors (default false)
% * 'pattern': struct with fields 'color', 'column_bd', 'rows'
%
%% read in optional arguments
defaults={'h',1e-5,'vectorized',false,'isdfx',false,'pattern',[]};
options=dde_set_options(defaults,varargin,'pass_on');
if ~options.isdfx
    dfx=@(x,v)dfxloc(f,x,v,options.h,options.vectorized);
else
    dfx=f;
end
pat=options.pattern;
%% if pattern not provided assume full pattern
nx=size(x,1);
if isempty(pat)
    nf=size(dfx(x,0*x),1);
    pat=dde_jacobianpattern(nx,sparse(ones(nf,nx)));
end
%% which deviations to take
nlindev=max(pat.deviation(:));
[idev1,idev2]=find(sparse(ones(nlindev)));
npairs=length(idev1);
lindev=full(min(1,sparse(pat.i_sorted(:,2),pat.deviation(:,1),...
    ones(size(pat.deviation,1),1))));
sel=idev1<=idev2;
dev1=lindev(:,idev1(sel));
dev2=lindev(:,idev2(sel));
nsel=sum(sel);
%% evaluate derivatives
if ~options.isdfx
    fdiffs=0.25*(dfx(x,dev1+dev2)-dfx(x,dev1-dev2));
else
    if options.vectorized
        xx=x(:,ones(1,nsel));
        fdiffs=0.25*(dfx(xx,dev1+dev2)-dfx(xx,dev1-dev2));
    else
        for i=ndev:-1:1
            fdiffs1(:,i)=dfx(x,dev1(:,i)+dev2(:,i));
            fdiffs2(:,i)=dfx(x,dev1(:,i)-dev2(:,i));
        end        
        fdiffs=0.25*(fdiffs1-fdiffs2);
    end
end
[~,lmap]=ismember([idev2(~sel),idev1(~sel)],[idev1,idev2],'rows');
fdiffall(:,sel)=fdiffs;
fdiffall(:,~sel)=fdiffall(:,lmap);
nf=size(fdiffall,1);
[~,lin_indev]=ismember(pat.deviation,[idev1,idev2],'rows');
%% insert values into sparse Jacobian
df.ind=pat.i_sorted;
lin_ind=sub2ind(size(fdiffall),pat.i_sorted(:,1),lin_indev);
df.vals=fdiffall(lin_ind);
%% remove exact zero entries
isnz=df.vals~=0;
df.vals=df.vals(isnz);
df.ind=df.ind(isnz,:);
%% duplicate symmetric entries
idiff=df.ind(:,2)<df.ind(:,3);
df.ind=[df.ind;df.ind(idiff,[1,3,2])];
df.vals=[df.vals;df.vals(idiff)];
[df.ind,ix]=sortrows(df.ind);
df.vals=df.vals(ix);
df.size=[nf,nx,nx];
end
%%
function fdiffs=dfxloc(f,x,dev,h,vectorized)
ndev=size(dev,2);
dev=dev*h;
xx=x(:,ones(1,ndev));
xu=xx+dev;
xl=xx-dev;
if vectorized
    fvals=f([xu,xl,x]);
    fu=fvals(:,1:ndev);
    fl=fvals(:,ndev+(1:ndev));
    fx=fvals(:,2*ndev+1);
else
    for i=ndev:-1:1
        fu(:,i)=f(xu(:,i));
        fl(:,i)=f(xl(:,i));
    end
    fx=f(x);
end
fdiffs=(fu+fl-2*fx(:,ones(1,ndev)))/h^2;
end
