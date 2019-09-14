function psol=dde_psol_from_x(x,template,free_par_ind)
%% Convert solution vector(s) into array of psol
%
% $Id$
%%
nps=size(x,2);
np=numel(template.profile);
psol=repmat(template,1,nps);
for i=nps:-1:1
    psol(i).profile(:)=x(1:np,i);
    psol(i).period=x(np+1,i);
    psol(i).parameter(free_par_ind)=x(np+2:end,i);
end
end
