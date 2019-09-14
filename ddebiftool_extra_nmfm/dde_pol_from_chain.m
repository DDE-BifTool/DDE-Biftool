function cfpol=dde_pol_from_chain(cf)
orders=unique(cf.gx);
nords=accumarray(cf.gx',ones(1,length(cf.gx)));
ind=length(cf.gx)+1;
for i=length(orders):-1:1
    cfpol(i).fac=[];
    cfpol(i).gx=orders(i);
    cfpol(i).mats=[];
    for j=1:nords(i)
        ind=ind-1;
        inz=find(cf.xo(:,ind));
        groups=cell(1,i);
        ic=0;
        for k=1:length(inz)
            groups{k}=ic+(1:cf.xo(inz(k),ind));
            ic=ic+cf.xo(inz(k),ind);
        end
        [mats,factors]=dde_polarization_coeffs(i,groups(1:k));
        matext=zeros(size(cf.xo,1),size(mats,2));
        matext(inz,:)=mats;
        cfpol(i).fac=[cfpol(i).fac,factors*cf.fac(ind)];
        cfpol(i).mats=[cfpol(i).mats,matext];
    end
    [mat,gc]=nmfm_matscale(cfpol(i).mats);
    [cfpol(i).mats,cfpol(i).fac]=nmfm_matcollect(mat,cfpol(i).fac.*gc);
end
end



