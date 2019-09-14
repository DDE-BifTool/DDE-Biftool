%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: tutorial_IV.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
clear variables
close all
addpath(...
    '../../../ddebiftool',...
    '../../../ddebiftool_utilities',...
    '../../../ddebiftool_extra_nmfm');

g=@(x) (tanh(x-1)+tanh(1))*cosh(1)^2;
neuron_sys_rhs = @(xx,par) [
        -xx(1,1,:)-par(1)*g(par(2)*xx(1,2,:))+par(3)*g(par(4)*xx(2,3,:));
        -xx(2,1,:)-par(1)*g(par(2)*xx(2,2,:))+par(3)*g(par(4)*xx(1,3,:))];
    
tau_ind=[5 6];
funcs=set_funcs(...
    'sys_rhs', neuron_sys_rhs,...
    'sys_tau', @() tau_ind,...
    'sys_deri', @neural_deri,...'sys_mfderi', @neural_mfderi,...
    'x_vectorized',true...
);

%% continue steady state
x0=[0; 0];

% use variable names instead of integers for parameter position
inda=1;
indc=3;

% set parameters
a=0.069;
b=2;
c=0.4;
d=1.2;
tau1=11.6;
tau2=20.3;

par=[a b c d tau1 tau2];

contpar=indc;
stst_branch=SetupStst(funcs,'x',x0,'parameter',par,'step',0.001,...
    'contpar',contpar,'max_step',[contpar,0.001],...
    'max_bound',[contpar,0.9]);

disp('--- Trivial equilibria ---');
figure(1);clf
stst_branch=br_contn(funcs,stst_branch,1000);
stst_branch=br_stabl(funcs,stst_branch,0,1);
nunst=GetStability(stst_branch);

% Bifurcation detection
[stst_branch,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch);

%% Plot bifurcation diagram
figure(2); clf
br_splot2(stst_branch,nunst,indc,'')
axis([0.4 0.9 -1 1])

%% Constructing an initial small-amplitude orbit near a Hopf bifurcation
flags_hopf=br_getflags(stst_branch,'hopf');
first_hopf=p_tohopf(funcs,stst_branch.point(flags_hopf(1)));
intervals=20;
degree=4;
[psol,stepcond]=p_topsol(funcs,first_hopf,1e-3,degree,intervals);
% correct periodic solution guess:
method=df_mthod(funcs,'psol');
[psol,success]=p_correc(funcs,psol,indc,stepcond,method.point)

%% add and plot stability
psol.stability=p_stabil(funcs,psol,method.stability);
figure(3);
p_splot(psol);

%% close-up of near 1
figure(4);
p_splot(psol);
axis([0.9999 1.0001   -0.0002    0.0002]);

%% manual construction and continuation of psol branch
per_orb=df_brnch(funcs,indc,'psol'); % empty branch:
deg_psol=p_topsol(funcs,first_hopf,0,degree,intervals);
per_orb.point=deg_psol;
%psol=rmfield(psol,'stability');
per_orb.point(2)=psol;
figure(5)
% compute periodic solutions branch
[per_orb,s,f,r]=br_contn(funcs,per_orb,53);
xlabel('c');ylabel('amplitude');

%% Branch off at first Hopf bifurcation
figure(6)
disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');
[per_orb,suc]=SetupPsol(funcs,stst_branch,flags_hopf(1),...
    'print_residual_info',0,'intervals',20,'degree',4,...
    'max_bound',[contpar,1],'max_step',[contpar,0.001]);
if ~suc
    error('initialization of periodic orbit failed');
end
per_orb=br_contn(funcs,per_orb,558);

% continue with smaller stepsize to detect PD4
per_orb.parameter.max_step=[3 1.0e-04];
per_orb=br_contn(funcs,per_orb,20);

% calculate stability
per_orb=br_stabl(funcs,per_orb,0,1);
[nunst_per,floqper,triv_defect,persols]=...
    GetStability(per_orb,'exclude_trivial',true);

%% One-parameter bifurcation diagram with the amplitude on the ordinate
figure(7)
br_splot2(stst_branch,nunst,indc,'');
hline = findobj(gcf, 'type', 'line');
set(hline,'Color',[.8 .8 .8])

br_splot2(per_orb,nunst_per,indc,'amplitude');
axis auto

text(0.702234157149676,0.423710735872078,'$NS_1$','Interpreter','latex')
text(0.643609307933753, 0.649755377784698,'$PD_1$','Interpreter','latex')
text(0.451, 2.00030155454736,'$LPC_1$','Interpreter','latex')
text(0.46, 2.22671664888698,'$PD_2$','Interpreter','latex')
text(0.59137261484202, 3.05293313349473,'$PD_3$','Interpreter','latex')
text(0.61736355481937, 2.97746143538152,'$LPC_2$','Interpreter','latex')
text(0.524634115649801, 0.196267005691675,'$PD_4$','Interpreter','latex')

%% One-parameter bifurcation diagram with max x_1(t) on the ordinate
figure(9)
br_splot2(stst_branch,nunst,indc,'')
hline = findobj(gcf, 'type', 'line');
set(hline,'Color',[.8 .8 .8])

% add labels of bifurcation points manually
text(0.7020,0.3178,'$NS_1$','Interpreter','latex')
text(0.6420,0.47,'$PD_1$','Interpreter','latex')
text(0.43,1.5204,'$LPC_1$','Interpreter','latex')
text(0.4420,1.77,'$PD_2$','Interpreter','latex')
text(0.575,2.58,'$PD_3$','Interpreter','latex')
text(0.6135,2.5382,'$LPC_2$','Interpreter','latex')
text(0.5220,1.4031,'$PD_4$','Interpreter','latex')

br_splot2(per_orb,nunst_per,indc,'max_x');
axis auto

%% periodic orbit to hopf point
hopf2=p_tohopf(funcs, per_orb.point(end));
[hopf2,suc]=p_correc(funcs,hopf2,indc,[],stst_branch.method.point)
hopf2.x
hopf2.stability=p_stabil(funcs,hopf2,stst_branch.method.stability);
hopf2.stability.l1(1:5)
hopf2=nmfm_hopf(funcs,hopf2);
fprintf('First Lyapunov coefficient l1: %g\n',hopf2.nmfm.L1);

%% Location of bifurcations points
fprintf('----- Location of bication points points ------ \n');

% Neimark-Sacker bifurcation
ind_ns=find(abs(diff(nunst_per))==2);
NS1=per_orb.point(ind_ns);
fprintf('Neimark-Sacker bifurcation detected near %f\n',...
    NS1.parameter(indc));

% fold and period doubling bifurcations
ind_bif=find(abs(diff(nunst_per))==1);
for i=1:length(ind_bif)
    p1=per_orb.point(ind_bif(i));
    
    % remove trivil multiplier
    [i1,i2]=min(abs(p1.stability.mu-1));
    p1.stability.mu(i2) = 0;
    [i1,i2]=min(abs(abs(p1.stability.mu)-1));
    s=sign(p1.stability.mu(i2));
    if s>0
        fprintf('Fold detected near %f\n',p1.parameter(indc));
    else
        fprintf('Period doubling detected near %f\n',p1.parameter(indc));
    end
end

%% bifurcation points (used in tutorial V)
PD1  = per_orb.point(ind_bif(1));
LPC1 = per_orb.point(ind_bif(2));
PD2  = per_orb.point(ind_bif(3));
PD3  = per_orb.point(ind_bif(4));
LPC2 = per_orb.point(ind_bif(5));
PD4  = per_orb.point(ind_bif(6));

%% Branch off at period doubling
figure(10)
[pd1,suc]=DoublePsol(funcs,per_orb,ind_bif(1));
if ~suc
    error('branching off at period doubling failed');
end
pd1.parameter.max_step=[indc 0.001];
pd1=br_contn(funcs,pd1,489);

pd1=br_stabl(funcs,pd1,0,0);
[nunst_per1,dom,triv_defect]=GetStability(pd1,'exclude_trivial',true);

%% One-parameter bifurcation diagram for family of periodic orbits
figure(11)
br_splot2(stst_branch,nunst,indc,'')
br_splot2(per_orb,nunst_per,indc,'max_x');
hline = findobj(gcf, 'type', 'line');
set(hline,'Color',[.8 .8 .8]);

br_splot2(pd1,nunst_per1,indc,'max_x',1.0);

figure(12)
br_splot2(stst_branch,nunst,indc,'')
br_splot2(per_orb,nunst_per,indc,'amplitude');
hline = findobj(gcf, 'type', 'line');
set(hline,'Color',[.8 .8 .8]);

br_splot2(pd1,nunst_per1,indc,'amplitude')

%% convert last point Hopf
hopf3=p_tohopf(funcs,pd1.point(end));
[hopf3,s]=p_correc(funcs,hopf3,indc,[],stst_branch.method.point)
hopf3.x
hopf3.stability=p_stabil(funcs,hopf3,stst_branch.method.stability);
hopf3.stability.l1(1:5)
% figure(4); clf
% p_splot(hopf3)
hopf3=nmfm_hopf(funcs,hopf3);
fprintf('First Lyapunov coefficient l1: %g\n',hopf3.nmfm.L1);

%% compare orbit on per_orb near PD_1 and first orbit on pd1
figure(13)
p_pplot(pd1.point(1),1);
hold on;
p_pplot(per_orb.point(ind_bif(1)));
ylabel('amplitude');

%% Branch off at PD_2
figure(14)
[pd2,suc]=DoublePsol(funcs,per_orb,ind_bif(3));
if ~suc
    error('branching off at period doubling failed');
end
pd2.parameter.max_step=[3 0.001];
pd2=br_contn(funcs,pd2,425);

% continue with smaller stepsize
pd2.parameter.max_step=[3 0.0001];
pd2=br_contn(funcs,pd2,33);

pd2=br_stabl(funcs,pd2,0,0);
[nunst_per2,dom,triv_defect]=GetStability(pd2,'exclude_trivial',true);

%% One-parameter bifurcation diagram for family of periodic orbits
% with ordinate max_x
figure(15)
br_splot2(stst_branch,nunst,indc,'')
br_splot2(per_orb,nunst_per,indc,'max_x',1.0);
hline = findobj(gcf, 'type', 'line');
set(hline,'Color',[.8 .8 .8])

br_splot2(pd2,nunst_per2,indc,'max_x')

% with ordinate amplitude
figure(16); clf
br_splot2(stst_branch,nunst,indc,'')
br_splot2(per_orb,nunst_per,indc,'amplitude',1.0);
hline = findobj(gcf, 'type', 'line');
set(hline,'Color',[.8 .8 .8])

br_splot2(pd2,nunst_per2,indc,'amplitude')

%% Branch off at PD3
figure(17)
[pd3,suc]=DoublePsol(funcs,per_orb,ind_bif(4),'radius',0.04);
if ~suc
    error('branching off at period doubling failed');
end

pd3.parameter.max_step=[3 0.001];

pd3=br_contn(funcs,pd3,230);
pd3=br_stabl(funcs,pd3,0,1);
[nunst_per3,dom,triv_defect]=GetStability(pd3,'exclude_trivial',true);

%% One-parameter bifurcation diagram for family of periodic orbits
% with ordinate amplitude
figure(18); clf
br_splot2(stst_branch,nunst,indc,'')
br_splot2(per_orb,nunst_per,indc,'amplitude',0.8);
hline = findobj(gcf, 'type', 'line');
set(hline,'Color',[.8 .8 .8])

br_splot2(pd3,nunst_per3,indc,'amplitude',0.8);

%% zoom in
figure(19); clf
br_splot2(stst_branch,nunst,indc,'')
br_splot2(per_orb,nunst_per,indc,'amplitude',0.8);
hline = findobj(gcf, 'type', 'line');
set(hline,'Color',[.8 .8 .8])
br_splot2(pd3,nunst_per3,indc,'amplitude',0.8);
axis([0.5567    0.6040    2.7668    3.0776]);
fig2=get(gcf,'Children');
set(fig2(1),'Location','NorthWest')

%% save workspace
remove_graphics()
%save('neural.mat')
