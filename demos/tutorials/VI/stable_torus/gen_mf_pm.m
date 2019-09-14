%% This file generates the multilinear forms needed for calculating 
%  the parameter-dependent normal form coefficients
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: gen_mf_pm.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
%% Define system
clear
n = 1; % number of components
taus = 9; % number of delays
free_pars = [1 2];

xx = sym('xx', [n,taus+1]);
par = sym('par', [1,13]);
taus = [0 par(5) par(6) par(7) par(8) par(9) par(10) par(11) par(12) par(13)];
m = length(taus);

% system
f = -par(4)*xx(1,1)-par(1)*xx(1,2)-par(2)*xx(1,3)-par(1)*par(3)*xx(1,1)*(par(1)*xx(1,4)+par(2)*xx(1,5)+par(4)*xx(1,2))-par(2)*par(3)*xx(1,1)*(par(1)*xx(1,5)+par(2)*xx(1,6)+par(4)*xx(1,3))-par(1)^2*par(3)^2*xx(1,1)*xx(1,2)*(par(1)*xx(1,7)+par(2)*xx(1,8)+par(4)*xx(1,4))-par(1)*par(2)*par(3)^2*xx(1,1)*xx(1,2)*(par(1)*xx(1,8)+par(2)*xx(1,9)+par(4)*xx(1,5))-par(2)*par(1)*par(3)^2*xx(1,1)*xx(1,3)*(par(1)*xx(1,8)+par(2)*xx(1,9)+par(4)*xx(1,5))-par(2)^2*par(3)^2*xx(1,1)*xx(1,3)*(par(1)*xx(1,9)+par(2)*xx(1,10)+par(4)*xx(1,6))-1/2*par(3)^2*xx(1,1)^2*par(1)*(par(4)^2*xx(1,2)+2*par(4)*(par(1)*xx(1,4)+par(2)*xx(1,5))+par(1)^2*xx(1,7)+2*par(1)*par(2)*xx(1,8)+par(2)^2*xx(1,9))-1/2*par(3)^2*xx(1,1)^2*par(2)*(par(4)^2*xx(1,3)+2*par(4)*(par(1)*xx(1,5)+par(2)*xx(1,6))+par(1)^2*xx(1,8)+2*par(1)*par(2)*xx(1,9)+par(2)^2*xx(1,10));

sys_dir='/home/maikel/Dropbox/Documents/Uni/Master/2015/Reading course DDE 2/VI/stable_torus/';

%% No changes needed from here. Just run the code.
%% generate free parameter file
fid = fopen('get_free_pars.m', 'w');
fprintf(fid, 'function free_pars = get_free_pars()\r\n');
fprintf(fid, 'free_pars=%s\r\n',mat2str(free_pars));
fprintf(fid, 'end');
fclose(fid);

%% B
phi0 = sym('phi0', [n*m,1]);
phi1 = sym('phi1', [n*m,1]);

Bphi1phi1=jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1;

matlabFunction(Bphi1phi1, 'file',strcat(sys_dir,'sys_B'),...
    'vars',{xx,par,phi0,phi1});

%% B2
phi0 = sym('phi0', [n*m,1]);
phi1 = sym('phi1', [n*m,1]);
psi0=sym('psi0',[2,1]);

Bphi1phi1psi0=jacobian(jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1,...
    [par(free_pars(1)); par(free_pars(2))])*psi0;

matlabFunction(Bphi1phi1psi0, 'file',strcat(sys_dir,'sys_B2'),...
    'vars',{xx,par,phi0,phi1,psi0});

%% C
phi0 = sym('phi0', [n*m,1]);
phi1 = sym('phi1', [n*m,1]);
phi2 = sym('phi2', [n*m,1]);
Bphi1phi1phi2=jacobian(jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1,xx(:))*phi2;
matlabFunction(Bphi1phi1phi2, 'file',strcat(sys_dir,'sys_C'),...
    'vars',{xx,par,phi0,phi1,phi2});

%% J1
J1=jacobian(f,[par(free_pars(1)); par(free_pars(2))]);
matlabFunction(J1, 'file',strcat(sys_dir,'sys_J1'), 'vars',{xx,par});

%% J2
psi0=sym('psi0',[2,1]);
psi1=sym('psi1',[2,1]);
J2=jacobian(jacobian(f,[par(free_pars(1)); par(free_pars(2))])*psi0,...
    [par(free_pars(1)); par(free_pars(2))])*psi1;
matlabFunction(J2, 'file',strcat(sys_dir,'sys_J2'), 'vars',{xx,par,psi0,psi1});

%% A1
phi0=sym('phi0', [n*m,1]);
psi1=sym('psi0',[2,1]);
A1=jacobian(jacobian(f,xx(:))*phi0,[par(free_pars(1)); par(free_pars(2))])*psi1;
fh = matlabFunction(A1, 'file',strcat(sys_dir,'sys_A1'),...
    'vars',{xx,par,phi0,psi0});

%% Big system
q0 = sym('q0',[n,1]);
q1 = sym('q1',[n,1]);
phi0 = sym('phi0', [n,m]);
phi1 = sym('phi1', [n,m]);
p1 = sym('p1',[n,1]);
assume(p1,'real');
p0 = sym('p0',[n,1]);
assume(p0,'real');
xi1 = sym('xi1',[n,1]);
xi2 = sym('xi2',[n,1]);
a2 = sym('a2');
b2 = sym('b2');
gamma1 = sym('gamma1');

% charisterictic matrix at lambda=0
D=zeros(n);
for j=1:m
    D=D-jacobian(f,xx(:,j));
end

p1sumBphi0=zeros(1,n);
for j=1:m
   p1sumBphi0=p1sumBphi0+p1'*jacobian(jacobian(f,xx(:,j))*phi0(:,j),xx(:,j));
end

p0sumBphi0=zeros(1,n);
for j=1:m
   p0sumBphi0=p0sumBphi0+p0'*jacobian(jacobian(f,xx(:,j))*phi0(:,j),xx(:,j));
end

p1sumBphi1=zeros(1,n);
for j=1:m
   p1sumBphi1=p1sumBphi1+p1'*jacobian(jacobian(f,xx(:,j))*phi1(:,j),xx(:,j));
end

J1=jacobian(f,[par(free_pars(1)); par(free_pars(2))]);

p1A1phi0=p1'*jacobian(jacobian(f,xx(:))*phi0(:),[par(free_pars(1)); par(free_pars(2))]);
p0A1phi0=p0'*jacobian(jacobian(f,xx(:))*phi0(:),[par(free_pars(1)); par(free_pars(2))]);
p1A1phi1=p1'*jacobian(jacobian(f,xx(:))*phi1(:),[par(free_pars(1)); par(free_pars(2))]);

bigsys=[-D, J1; p1sumBphi0, p1A1phi0; p0sumBphi0+p1sumBphi1, p0A1phi0+p1A1phi1];


matlabFunction(bigsys, 'file',strcat(sys_dir,'sys_bigsys'),...
    'vars',{xx,par,p1,p0,phi1,phi0});

% rhs bigsys
% derivatives charisterictic matrix at lambda=0
dD=eye(n);
ddD=zeros(n);
dddD=zeros(n);
ddddD=zeros(n);
dddddD=zeros(n);
for j=1:m
    dD=dD+taus(j)*jacobian(f,xx(:,j));
    ddD=ddD-taus(j)^2*jacobian(f,xx(:,j));
    dddD=dddD+taus(j)^3*jacobian(f,xx(:,j));
    ddddD=ddddD-taus(j)^4*jacobian(f,xx(:,j));
    dddddD=dddddD+taus(j)^5*jacobian(f,xx(:,j));
end
L11=dD*q1+1/2*ddD*q0;

Bphi1phi1=jacobian(jacobian(f,xx(:))*phi1(:),xx(:))*phi1(:);
p1sumBphi0H0010=0;
for j=2:m
   p1sumBphi0H0010=p1sumBphi0H0010+p1'*(jacobian(jacobian(...
       f,xx(:,j))*phi0(:,j),xx(:,j))*(taus(j)*q1+1/2*taus(j)*q0));
end

L21=1/2*p1'*Bphi1phi1+p1sumBphi0H0010;

phi0sunH1100=p0'*( (b2/6*dddD+a2/12*ddddD+ddD/2*gamma1)*q0 ...
                + (b2/2*ddD+a2/3*dddD)*q1 + ddD/2*xi1 + dD*xi2) ...
             + p1'*( (b2/24*ddddD-a2/60*dddddD+gamma1/6*dddD)*q0 ...
                + (b2/6*dddD+a2/12*ddddD)*q1 + dddD/6*xi1 + ddD/2*xi2);

p0sumBphi0H0010=0;
p1sumBphi1H0010=0;
for j=2:m
   p0sumBphi0H0010=p0sumBphi0H0010+p0'*(jacobian(jacobian(...
       f,xx(:,j))*phi0(:,j),xx(:,j))*(taus(j)*q1+1/2*taus(j)^2*q0));
   p1sumBphi1H0010=p1sumBphi1H0010+p1'*(jacobian(jacobian(...
       f,xx(:,j))*phi1(:,j),xx(:,j))*(taus(j)*q1+1/2*taus(j)^2*q0));
end   

L31=3*phi0sunH1100-p0'*Bphi1phi1-p0sumBphi0H0010-p1sumBphi1H0010;
            
bigsysRHS=[L11, zeros(2,1); L21, 0; L31, 1];

matlabFunction(bigsysRHS, 'file',strcat(sys_dir,'sys_bigsysRHS'),...
    'vars',{xx,par,q0,q1,p1,p0,phi1,phi0,a2,b2,gamma1,xi1,xi2});
