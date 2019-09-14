%% pseudo-vectorize/reshape sys_dirderi
%
% $Id: dde_wrap_dirderi.m 348 2019-06-19 13:09:01Z jansieber $
%%
function J=dde_wrap_dirderi(x,p,dx,dp,dirderi_arg,isxvec,ispvec,itau)
xsize=size(x);
psize=size(p);
vecdim=[xsize(3:end),1];
x=reshape(x,xsize(1),xsize(2),[]);
p=reshape(p,psize(1),psize(2),[]);
dx=reshape(dx,size(x));
dp=reshape(dp,psize(1),psize(2),[]);
npvec=size(p,3);
nvec=size(x,3);
if nargin<=7
    dirderi=dirderi_arg;
    Jrowdim=xsize(1);
else
    dirderi=@(xa,pa,dxa,dpa)dirderi_arg(itau,xa,pa,dxa,dpa);
    Jrowdim=1;
end
if nvec==1 || (isxvec && ~isxvec &&  npvec==1) || (isxvec && ispvec)
    %% complete vectorization
    J=dirderi(x,p,dx,dp);
elseif  (npvec==1 ||(all(all(diff(dp,[],3)==0)) && all(all(diff(p,[],3)==0)))) && (isxvec || nvec==1)
    %% in x vectorized, all p and dp are equal
    J=dirderi(x,p(1,:,1),dx,dp(1,:,1));
elseif isxvec
    %% in x vectorized, hopefully some sequential p and dp are equal
    pd=abs(diff(p,[],3))+abs(diff(dp,[],3));
    inddiff=find(any(pd,2));
    ibdlow=[1;inddiff+1];
    ibdup=[inddiff;size(p,3)];
    for i=length(ibdlow):-1:1
        rg=ibdlow(i):ibdup(i);
        ip=mod(rg(1)-1,npvec)+1;
        J(:,rg)=dirderi(x(:,:,rg),p(1,:,ip),dx(:,:,rg),dp(1,:,ip));
    end
else
    for i=size(x,3):-1:1
        ip=mod(i-1,npvec)+1;
        J(:,i)=dirderi(x(:,:,i),p(1,:,ip),dx(:,:,i),dp(1,:,ip));
    end
end
J=reshape(J,[Jrowdim,vecdim]);
end
