function x=dde_x_from_hopf(hopf,free_par_ind)
%% Convert hopf structure array into vector(s)
%
% $Id$
%%
if isempty(hopf)
    x=[];
    return
end
par=NaN(length(free_par_ind),length(hopf));
for i=length(hopf):-1:1
    xval(:,i)=hopf(i).x;
    v(:,i)=hopf(i).v;
    par(:,i)=hopf(i).parameter(free_par_ind).';
    omega(i)=hopf(i).omega;
end
x=[xval;real(v);imag(v);omega;par];
end