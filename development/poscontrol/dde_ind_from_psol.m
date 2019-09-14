function ind=dde_ind_from_psol(psol,free_par_ind)
%% Determine indices for entries in vector from structure
%
% $Id$
%%
ind=psol(1);
ind.profile=reshape(1:numel(psol.profile),size(psol.profile));
ind.period=ind.profile(end)+1;
ind.parameter=ind.period+(1:length(free_par_ind));
end