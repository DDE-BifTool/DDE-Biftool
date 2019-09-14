function [w,base_v,Dout]=coll_barywt(mesh,degree,difforder)
npi=degree+1;
base_h=mesh(1:npi)/mesh(npi)*2-1;
base_v=base_h';
xdiff=base_h(ones(npi,1),:)-base_v(:,ones(npi,1));
w=1./prod(eye(npi)-xdiff,2);
Dout=eye(npi);
if difforder>0
    wrep=w(:,ones(npi,1));
    denom=(xdiff'-xdiff);
    denom(1:npi+1:end)=Inf;
    D=wrep'./wrep./denom;
    Drsum=sum(D,2);
    D(1:npi+1:end)=-Drsum;
    for i=1:difforder
        Dout=4*D*Dout;
    end
end
end
