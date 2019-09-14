function [iout,len]=dde_ind_from_psolbif(pt,free_par,base,vfnames)
[psol,cv]=dde_coll_convert(pt);
ifree=[free_par,cv.extra_ind.parameter];
[indpsol,len]=dde_ind_from_point(psol,ifree,base);
indpsol_par=indpsol.parameter;
indpsol.parameter=zeros(1,length(indpsol.parameter));
indpsol.parameter(ifree)=indpsol_par;
ind=dde_coll_convert(indpsol,cv);
ind.parameter=ind.parameter(ind.parameter>0);
for i=1:length(vfnames)
    iout.(vfnames{i})=ind.(vfnames{i});
end
end
