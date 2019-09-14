function psol=psol_from_x(x,template,free_par_ind)
npoints=size(x,2);
psol=repmat(template(1),npoints,1);
np=numel(template(1).profile);
for i=1:npoints
    psol(i).profile(:)=x(1:np,i);
    psol(i).period=x(np+1,i);
    psol(i).parameter(free_par_ind)=x(np+2:end,i);
end
end