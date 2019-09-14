%% pseudo-vectorize/reshape sys_rhs
%
% $Id: dde_wrap_rhs.m 348 2019-06-19 13:09:01Z jansieber $
%%
function y=dde_wrap_rhs(x,p,rhs,isxvec,ispvec,itau)
xsize=size(x);
psize=size(p);
vecdim=[xsize(3:end),1];
x=reshape(x,xsize(1),xsize(2),[]);
p=reshape(p,psize(1),psize(2),[]);
npvec=size(p,3);
nvec=size(x,3);
if nargin<=5
    locrhs=rhs;
    yrowdim=xsize(1);
else
    locrhs=@(xa,pa)rhs(itau,xa,pa);
    yrowdim=1;
end
if nvec==1 || (isxvec && ~isxvec && npvec==1) || (isxvec && ispvec)
    y=locrhs(x,p);
elseif  (npvec==1 || all(all(diff(p,[],3)==0))) && (isxvec || nvec==1)
    y=locrhs(x,p(1,:,1));
elseif isxvec
    %% in x vectorized, hopefully some sequential p are equal
    pd=diff(p,[],3);
    inddiff=find(any(pd,2));
    ibdlow=[1;inddiff+1];
    ibdup=[inddiff;size(p,3)];
    for i=length(ibdlow):-1:1
        rg=ibdlow(i):ibdup(i);
        ip=mod(rg(1)-1,npvec)+1;
        y(:,rg)=locrhs(x(:,:,rg),p(1,:,ip));
    end
else
    for i=size(x,3):-1:1
        ip=mod(i-1,npvec)+1;
        y(:,i)=locrhs(x(:,:,i),p(1,:,ip));
    end
end
y=reshape(y,[yrowdim,vecdim]);
end
