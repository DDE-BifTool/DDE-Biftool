%% reshape deri from vectorization  or pseudo-vectorize
%
% $Id: dde_wrap_deri.m 348 2019-06-19 13:09:01Z jansieber $
%%
function J=dde_wrap_deri(x,p,nx,np,v,sys_deri_arg,isxvec,ispvec,itau)
xsize=size(x);
psize=size(p);
x=reshape(x,xsize(1),xsize(2),[]);
p=reshape(p,psize(1),psize(2),[]);
npvec=size(p,3);
nvec=size(x,3);
if isempty(v)
    v=zeros(0,size(x,3));
end
if nargin<=8
    sys_deri=sys_deri_arg;
    istau=false;
else
    sys_deri=@(xa,pa,nxa,npa,v)sys_deri_arg(itau,xa,pa,nxa,npa);
    istau=true;
end
if nvec==1 || (isxvec && ~isxvec && npvec==1) || (isxvec && ispvec)
    J=sys_deri(x,p,nx,np,v);
elseif (npvec==1 || all(all(diff(p,[],3)==0))) && (isxvec || nvec==1)
    J=sys_deri(x,p(1,:,1),nx,np,v);
elseif isxvec
    %% in x vectorized, hopefully some sequential p and dp are equal
    pd=abs(diff(p,[],3))+abs(diff(dp,[],3));
    inddiff=find(any(pd,2));
    ibdlow=[1;inddiff+1];
    ibdup=[inddiff;size(p,3)];
    for i=length(ibdlow):-1:1
        rg=ibdlow(i):ibdup(i);
        ip=mod(rg(1)-1,npvec)+1;
        J(:,:,rg)=sys_deri(x(:,:,rg),p(1,:,ip),nx,np,v(:,rg));
    end
else
    for i=nvec:-1:1
        ip=mod(i-1,npvec)+1;
        J(:,:,i)=sys_deri(x(:,:,i),p(1,:,ip),nx,np,v(:,i));
    end
end
if ~istau
    J=reshape(J,size(x,1),[],nvec);
else
   if length(nx)==1 && isempty(np)
       J=reshape(J,[1,xsize(1),nvec]);
   elseif length(nx)==2
       J=reshape(J,[xsize(1),xsize(1),nvec]);
   elseif  length(nx)==1 && length(np)==1
       J=reshape(J,[1,xsize(1),nvec]);
   end
end
end

