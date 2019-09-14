function bt=p_tobt(funcs,point)
%% convert point to Bogdanov-Takens bifurcation point
% function bt=p_tobt(funcs,point)
% INPUT:
%   funcs : problem functions
%	point : point approximating the BT-point
% OUTPUT:
%	bt : uncorrected starting guess for bt point
%
% $Id: p_tobt.m 309 2018-10-28 19:02:42Z jansieber $
%%
bt=point;

bt.kind='BT';
bt.parameter=point.parameter;

%% set up bordered system
n=length(bt.x);
D=ch_matrix(funcs,bt.x,bt.parameter,0);
dD=ch_matrix(funcs,bt.x,bt.parameter,0,'deri',1);

[X,L] = eig(D);
[~,I] = sort(diag(L));
c0=real(X(:,I(1)));
[X,L] = eig(D');
[~,I] = sort(diag(L));
b0=real(X(:,I(1)));

Bord = [ D b0; c0' 0];

%% calculate vectors q0 and p1
q0=Bord\[zeros(n,1);1];
q0=q0(1:n);
q0=q0/norm(q0);
q1 = Bord\[-dD*q0; 0];
q1 = q1(1:n);
q1=q1-(q1'*q0)*q0;
%% calculate transpose vectors p1 and p0
Bord=Bord';

p1=Bord\[zeros(n,1);1];
p1=p1(1:n);

p0 = Bord\[-dD'*p1; 0];
p0 = p0(1:n);

%% return the vectors
bt.nvec.p1=p1;
bt.nvec.p0=p0;
bt.nvec.q0=q0;
bt.nvec.q1=q1;
if isfield(bt,'omega')
    bt.omega=0;
end
if isfield(bt,'v')
    bt.v=q0;
end

end
