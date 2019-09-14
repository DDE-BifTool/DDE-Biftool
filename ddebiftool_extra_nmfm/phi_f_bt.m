function values = phi_f_bt(funcs, point)
%% detection function for Bogdanov-Takens bifurcation along fold curve
%
% $Id: phi_f_bt.m 309 2018-10-28 19:02:42Z jansieber $
%% 
np=length(point);
values=NaN(size(point));
D=ch_matrix(funcs,point(1).x,point(1).parameter,0);
[E1,E2]=eig(D');
[i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
p1=real(E1(:,i2));
p1old=p1;
for i=1:np
    D=ch_matrix(funcs,point(i).x,point(i).parameter,0);
    [E1,E2]=eig(D');
    [i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
    p1=real(E1(:,i2));
    p1=p1*sign(p1'*p1old);
    dD=ch_matrix(funcs,point(i).x,point(i).parameter,0,'deri',1);
    values(i)=p1'*dD*point(i).v;
end
end
