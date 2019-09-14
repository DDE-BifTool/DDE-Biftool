%% Stable torus in approxiamtion to DDESD with constant delays
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: stable_torus.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
clear
close all
load 'stable_torus.mat'
%% Define system
clear
close all
addpath(...
    '../../../../ddebiftool',...
    '../../../../ddebiftool_utilities',...
    '../../../../ddebiftool_extra_psol',...
    '../../../../ddebiftool_extra_nmfm',...
    '../../IV');

sys_rhs = @(xx,par) -par(4).*xx(1,1,:)-par(1).*xx(1,2,:)-par(2).*xx(1,3,:)...
    -par(1).*par(3).*xx(1,1,:).*(par(1).*xx(1,4,:)...
    +par(2).*xx(1,5,:)+par(4).*xx(1,2,:))...
    -par(2).*par(3).*xx(1,1,:).*(par(1).*xx(1,5,:)...
    +par(2).*xx(1,6,:)+par(4).*xx(1,3,:))...
    -par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,7,:)...
    +par(2).*xx(1,8,:)+par(4).*xx(1,4,:))...
    -par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,8,:)...
    +par(2).*xx(1,9,:)+par(4).*xx(1,5,:))...
    -par(2).*par(1).*par(3).^2.*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,8,:)...
    +par(2).*xx(1,9,:)+par(4).*xx(1,5,:))...
    -par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,9,:)...
    +par(2).*xx(1,10,:)+par(4).*xx(1,6,:))...
    -1/2.*par(3).^2.*xx(1,1,:).^2.*par(1).*(par(4).^2.*xx(1,2,:)...
    +2.*par(4).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:))...
    +par(1).^2.*xx(1,7,:)+2.*par(1).*par(2).*xx(1,8,:)...
    +par(2).^2.*xx(1,9,:))...
    -1/2.*par(3).^2.*xx(1,1,:).^2.*par(2).*(par(4).^2.*xx(1,3,:)...
    +2.*par(4).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:))...
    +par(1).^2.*xx(1,8,:)+2.*par(1).*par(2).*xx(1,9,:)...
    +par(2).^2.*xx(1,10,:));


tau_ind=[5 6 7 8 9 10 11 12 13];
funcs=set_funcs(...
    'sys_rhs', sys_rhs,...
    'sys_tau', @() tau_ind,...
    'sys_deri',@(xx,par,nx,np,v) sys_deri(xx,par,nx,np,v),...
    'sys_mfderi',@(xx,par,varargin) sys_mfderi(xx,par,varargin{:}),...
    'x_vectorized',1);

%% steady state
x0=0;

% use variable names instead of integers for parameter position
ind_kappa1=1;
ind_kappa2=2;

% set parameters
kappa1=2.6;
kappa2=3.4;
c=1;
gamma=4.75;
a1=1.3;
a2=6;

tau1=a1;
tau2=a2;
tau3=2*a1;
tau4=a1+a2;
tau5=2*a2;
tau6=3*a1;
tau7=a2+2*a1;
tau8=2*a2+a1;
tau9=3*a2;

par=[kappa1 kappa2 c gamma tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9];

stst.kind='stst';
stst.parameter=par;
stst.x=x0;

method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability)
figure(1);clf
p_splot(stst)
%% steady-state branch
kappa1=2.6;
kappa2=3;

contpar=ind_kappa2;
% setup steady-state branch
stst_br=SetupStst(funcs,'x',x0,'parameter',par,'step',0.01,...
    'contpar',contpar,'max_step',[contpar,0.01],...
    'max_bound',[[ind_kappa1,3.2]; [ind_kappa2, 4]],...
    'min_bound',[[ind_kappa1,1.6]; [ind_kappa2, 3]]);

disp('--- Trivial equilibria ---');
figure(1);clf
stst_br=br_contn(funcs,stst_br,1000);

% stability calculation
stst_br=br_stabl(funcs,stst_br,0,0);
nunst=GetStability(stst_br);

% bifurcation detection
[stst_br,stst_testfuncs]=LocateSpecialPoints(funcs,stst_br);

%% continue Hopf points in (kappa1,kapp2)
hopf_ind=br_getflags(stst_br,'hopf');
fprintf('----- Hopf branch 1 -----\n');
[hopf_branch1,suc] = SetupHopf(funcs, stst_br, hopf_ind(1),...
    'contpar', [ind_kappa1,ind_kappa2],...
    'dir', ind_kappa1);

hopf_branch1.parameter.max_step=[];

figure(2);clf; hold on
hopf_branch1=br_contn(funcs,hopf_branch1,300);
hopf_branch1=br_rvers(hopf_branch1);
hopf_branch1=br_contn(funcs,hopf_branch1,300);

fprintf('----- Hopf branch 2 -----\n');
[hopf_branch2,suc] = SetupHopf(funcs, stst_br, hopf_ind(2),...
    'contpar', [ind_kappa1,ind_kappa2],...
    'dir', ind_kappa1);

hopf_branch2.parameter.max_step=[];

hopf_branch2=br_contn(funcs,hopf_branch2,300);
hopf_branch2=br_rvers(hopf_branch2);
hopf_branch2=br_contn(funcs,hopf_branch2,300);

% detect codim-2 points
% hopf_branch1.method.bifurcation.secant_tolerance=1.0e-32;
% hopf_branch2.method.bifurcation.secant_tolerance=1.0e-32;
% hopf_branch1.method.bifurcation.secant_iterations=100;
% hopf_branch2.method.bifurcation.secant_iterations=100;

hopf_branch1.method.point.newton_max_iterations=100;
hopf_branch1.method.point.newton_nmon_iterations=100;
hopf_branch1.method.point.minimal_accuracy=1.0e-14;
hopf_branch1.method.point.halting_accuracy=1.0e-14;

hopf_branch2.method.point.newton_max_iterations=100;
hopf_branch2.method.point.newton_nmon_iterations=100;
hopf_branch2.method.point.minimal_accuracy=1.0e-14;
hopf_branch2.method.point.halting_accuracy=1.0e-14;

fprintf('----- Codimension-two detection along the first Hopf branch -----\n');
[hopf_branch1,hopf1_testfuncs]=LocateSpecialPoints(funcs,hopf_branch1);

fprintf('----- Codimension-two detection along the second Hopf branch -----\n');
[hopf_branch2,hopf2_testfuncs]=LocateSpecialPoints(funcs,hopf_branch2);

hoho_ind=br_getflags(hopf_branch1,'hoho');
hopf_branch1.point(hoho_ind).parameter(1:2)

hoho_ind=br_getflags(hopf_branch2,'hoho');
hopf_branch2.point(hoho_ind).parameter(1:2)
%% plot
figure(2);clf; hold on
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br);
plot(getpars(hopf_branch1.point,ind_kappa1),getpars(hopf_branch1.point,ind_kappa2))
plot(getpars(hopf_branch2.point,ind_kappa1),getpars(hopf_branch2.point,ind_kappa2))
legend('hopf\_branch1','hopf\_branch2');
xlabel('$\kappa_1$','Interpreter','LaTex');
ylabel('$\kappa_2$','Interpreter','LaTex');

%% Branch off at hopf near Hopf-Hopf on the branch hopf_branch1
hoho_ind=br_getflags(hopf_branch1,'hoho');
disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');
[per_orb,suc]=SetupPsol(funcs,hopf_branch1,hoho_ind-20,...
    'intervals',20,'degree',4,'contpar',ind_kappa1);
if ~suc
    error('initialization of periodic orbit failed');
end
figure(3); clf;
per_orb.parameter.max_step=[[ind_kappa1 0.002]; [ind_kappa2 0.002]];
per_orb.parameter.free=ind_kappa2;
per_orb=br_contn(funcs,per_orb,100);
per_orb=br_stabl(funcs,per_orb,0,1);
[nunst_per]=GetStability(per_orb,'exclude_trivial',true);
% nunst_per'
clf;br_splot2(per_orb,nunst_per,ind_kappa2,'amplitude');
xlabel('$\kappa_2$','Interpreter','LaTex');
ylabel('amplitude');

%% Branch off at hopf near Hopf-Hopf on the branch hopf_branch2
hoho_ind=br_getflags(hopf_branch2,'hoho');
disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');
[per_orb2,suc]=SetupPsol(funcs,hopf_branch2,hoho_ind-16,...
    'intervals',20,'degree',4,'contpar',ind_kappa1);
if ~suc
    error('initialization of periodic orbit failed');
end
figure(4); clf;
hold on
per_orb2.parameter.free=ind_kappa2;
per_orb2.parameter.max_step=[ind_kappa2 0.002];
per_orb2=br_contn(funcs,per_orb2,7);
per_orb2=br_stabl(funcs,per_orb2,0,1);
[nunst_per2]=GetStability(per_orb2,'exclude_trivial',true);
% nunst_per2'
clf;br_splot2(per_orb2,nunst_per2,ind_kappa2,'amplitude');
xlabel('$\kappa_2$','Interpreter','LaTex');
ylabel('amplitude');

%% continue Neimark-Sacker points
%% continue NS1 point in (kappa1,kappa2)
figure(5); clf
ind_ns=find(abs(diff(nunst_per))==2);
[ns1_pfuncs,ns1_br]=SetupTorusBifurcation(funcs,per_orb,ind_ns,...
    'contpar',[ind_kappa1,ind_kappa2],'dir',ind_kappa1,'step',1e-2);
ns1_br.parameter.max_step=[];
ns1_br=br_contn(ns1_pfuncs,ns1_br,50);
ns1_br=br_rvers(ns1_br);
ns1_br=br_contn(ns1_pfuncs,ns1_br,22);

%% continue NS2 point in (kappa1,kappa2)
figure(6); clf
ind_ns2=find(abs(diff(nunst_per2))==2);
[ns2_pfuncs,ns2_br]=SetupTorusBifurcation(funcs,per_orb2,ind_ns2(1),...
    'contpar',[ind_kappa1,ind_kappa2],'dir',ind_kappa1,'step',1e-2);
ns2_br.parameter.max_step=[];
ns2_br=br_contn(ns2_pfuncs,ns2_br,65);
ns2_br=br_rvers(ns2_br);
ns2_br=br_contn(ns1_pfuncs,ns2_br,10);
%% plot
figure(7);clf;hold on
cm=colormap('lines');
hopf_branch1_pl=plot(getpars(hopf_branch1.point,ind_kappa1),getpars(hopf_branch1.point,ind_kappa2),'color',cm(1,:))
plot(getpars(hopf_branch2.point,ind_kappa1),getpars(hopf_branch2.point,ind_kappa2),'color',cm(1,:))
ns1_br_pl=plot(getpars(ns1_br.point,ind_kappa1),getpars(ns1_br.point,ind_kappa2),'color',cm(7,:))
plot(getpars(ns2_br.point,ind_kappa1),getpars(ns2_br.point,ind_kappa2),'color',cm(7,:))

legend([hopf_branch1_pl ns1_br_pl],{'Hopf branches','Neimark-Sacker branches'})
xlabel('$\kappa_1$','Interpreter','LaTex');
ylabel('$\kappa_2$','Interpreter','LaTex');

%% investigate near neimark-sacker cusp
ns3_br=ns2_br;
ns3_br.point=ns2_br.point(14:15); % points near the turning point
ns3_br.parameter.max_step=[ind_kappa1 0.00005];

figure(8); clf
ns3_br=br_contn(ns2_pfuncs,ns3_br,130);
%% compare ns2_br and ns3_br
figure(9); clf; hold on
plot(getpars(ns2_br.point,ind_kappa1),getpars(ns2_br.point,ind_kappa2),'color',cm(1,:))
plot(getpars(ns3_br.point,ind_kappa1),getpars(ns3_br.point,ind_kappa2),'color',cm(7,:))
axis([2.502    2.504    3.5984    3.5994])
legend('ns2\_br','ns3\_br')
xlabel('$\kappa_1$','Interpreter','LaTex');
ylabel('$\kappa_2$','Interpreter','LaTex');

%% save workspace
save 'stable_torus.mat'
