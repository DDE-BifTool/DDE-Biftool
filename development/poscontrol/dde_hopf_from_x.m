function hopf=dde_hopf_from_x(x,template,free_par_ind)
%% Convert solution vector(s) into array of hopf points
%
% $Id$
%%
nhopfs=size(x,2);
hopf=repmat(template,1,nhopfs);
if isempty(x)
    return
end
n=size(hopf(1).x,1);
for i=1:nhopfs
    hopf(i).x=x(1:n,i);
    hopf(i).v=x(n+(1:n),i)+1i*x(2*n+(1:n),i);
    hopf(i).omega=x(3*n+1,i);
    hopf(i).parameter(free_par_ind)=...
        x(3*n+1+(1:length(free_par_ind)),i);
end
end
