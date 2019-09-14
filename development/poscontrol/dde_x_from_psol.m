function x=dde_x_from_psol(psol,free_par_ind)
%% Convert psol structure into vector
%
% $Id$
%%
if isempty(psol)
    x=[];
    return
end
par=NaN(length(free_par_ind),length(psol));
for i=length(psol):-1:1
    par(:,i)=psol(i).parameter(free_par_ind).';
    T(i)=psol(i).period;
    traj(:,i)=psol(i).profile(:);
end
x=[traj;T(:).';par];
end