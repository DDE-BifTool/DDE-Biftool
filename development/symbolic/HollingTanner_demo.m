%% Holling-Tanner model - Demo for Takens-Bogdanov point normal form
%
%% Differential equations for predator-prey model
%
% From:
% Liu, X., Liu, Y., and Wang, J. (2013). Bogdanov-Takens bifurcation of a 
% delayed ratio-dependent Holling-Tanner predator prey system. 
% In Abstract and Applied Analysis, volume 2013
% 
% $$x'=(x+m)(1-x-m)-\frac{xy}{ay+x}$$ 
%
% $$y'=\delta y\left(\beta-\frac{y(t-\tau)}{x(t-\tau)}\right)$$
% 
%
% <html>
% $Id: HollingTanner_demo.m 174 2017-03-10 21:59:17Z jansieber $
% </html>
%%
%#ok<*UNRCH>
clear;                           % clear variables
close all;                       % close figures
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');
%% Right-hand side
% Change the value of use_symbolic_results to false to use hand-coded
% right-hand side and derivatives. The structure ind contains indices at
% which the corresponding parameter with name from parnames is in the
% parameter vector.
ntau=1;
parnames={'beta','tau','a','m','h','delta'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% This problem definition relies on the file sym_Holling.m 
% created by the symbolic toolbox
fsym=set_symfuncs(@sym_Holling,'sys_tau', @()ind.tau);
%% Manually coded right-hand sides
fmanual=set_funcs(...
    'sys_rhs', @HollingTanner_rhs,...
    'sys_tau', @()ind.tau,...
    'sys_deri', @HollingTanner_deri,...
    'sys_mfderi',@HollingTanner_mfderi);
%% Define only r.h.s and a single directional derivative
fdir=set_funcs(...
    'sys_rhs', @HollingTanner_rhs,...
    'sys_tau', @()ind.tau,...
    'sys_dirderi',{@(x,p,dx,dp)dir_Holling_rhs_1(x,p,dx,dp,ind)},'x_vectorized',true);
%% Choose which problem definition is tested
funcs=fsym;
%% Sequence of parameters
getpar=@(x,i)arrayfun(@(p)p.parameter(i),x.point);
getx=@(x,i)arrayfun(@(p)p.x(i),x.point);
bgetpar=@(x,i,bif)arrayfun(@(p)p.parameter(i),x.point(br_getflags(x,bif)));
bgetx=@(x,i,bif)arrayfun(@(p)p.x(i),x.point(br_getflags(x,bif)));
%% Continuation of steady state branch
fprintf('----- Steady-state branch -----\n');
beta=0.5;
a=0.5; 
m=(1/30)*(1-beta/(a*beta+1)); 
h=(1/4)*(beta/(a*beta+1)-1)^2+m*beta/(a*beta+1); 
tau=1/4*(a*beta+1)^2/beta; 
beta=0.4; 
delta=0.5409;
parameter=cellfun(@(x)evalin('caller',x),parnames); %[beta,tau,a,m,h,delta];
parameter_bd={'max_bound',[ind.beta,0.6; ind.delta,0.7],...
    'min_bound',[ind.beta,0.4;ind.delta,0.4],...
    'max_step',[0,0.1;ind.beta,5e-3;ind.delta,5e-3]};
xster=-(1/2)*((beta/(a*beta+1)+2*m-1)+sqrt((1-2*m-beta/(a*beta+1))^2+4*(m*(1-m)-h)));
yster=beta*xster;

contpar=ind.beta;
stst_branch0 = SetupStst(funcs,'x',[xster; yster],'parameter',parameter,...
    'contpar',contpar,'max_step',[0,0.005],'min_bound',...
    [contpar 0.4],'max_bound',[contpar 0.6],...
    'newheuristics_tests',0);
figure(1);clf
ax1=gca;
title(ax1,sprintf('steady states for delta=%g',parameter(ind.delta)));
[stst_branch0] = br_contn(funcs,stst_branch0,100);
[stst_branch_wbifs,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch0);
nunst_stst=GetStability(stst_branch_wbifs);
%% Plot bifurcation diagram
beta_stst=getpar(stst_branch_wbifs,ind.beta);
x1_stst=getx(stst_branch_wbifs,1);
cla(ax1);
Plot2dBranch(stst_branch_wbifs,'y',@(p)p.x(1),'ax',ax1);
xlabel(ax1,'\beta');
ylabel(ax1,'x_1');
title(ax1,sprintf('steady states for delta=%g',parameter(ind.delta)));
%% Continuation of Hopf bifurcation in two parameters
% after initialization of hopfbranch
fprintf('----- Hopf branch -----\n');
[hopf_branch0,suc] = SetupHopf(funcs, stst_branch_wbifs, br_getflags(stst_branch_wbifs,'hopf'),...
    'contpar', [ind.delta,ind.beta],...
    'dir', ind.beta, 'step', 0.002,parameter_bd{:});
disp(['suc=',num2str(suc)]);
figure(2);clf
ax2=gca;
title(ax2,'Hopf in beta-delta plane');
hopf_branch0=br_contn(funcs,hopf_branch0,300);
%% Detect special points along Hopf curve
fprintf('----- Codimension-two detection along Hopf branch -----\n');
[hopf_branch_wbifs,hopftestfuncs]=LocateSpecialPoints(funcs,hopf_branch0);
nunst_hopf=GetStability(hopf_branch_wbifs,'exclude_trivial',true);
%% Singularity of defining system for Hopf bifurcation 
% Let us demonstrate that the defining system for the Hopf bifurcation is
% singular in the Takens-Bogdanov bifurcation. The 8x9 matrix |J0| below is its
% Jacobian (without the pseudo-arclength condition). Its singular-value
% decomposition shows one zero (up to round-off). |hbt| is the point along
% the Hopf branch closest to the BT point.
%
% The characteristic matrix has a double root at lambda=0 in BT.
hbt=hopf_branch_wbifs.point(br_getflags(hopf_branch_wbifs,'BT'));
J0=hopf_jac(funcs,hbt.x,hbt.omega,hbt.v,hbt.parameter,[ind.delta,ind.beta],hbt.v.');
disp('J0=');
disp(J0);
fprintf('Size of J0: (%d x %d), min(svd(J0))=%g\n',...
    size(J0,1),size(J0,2),min(svd(J0)));
charfunBT=@(lambda)det(ch_matrix(funcs,hbt.x(:,[1,1]),hbt.parameter,lambda));
figure(5);clf;
fplot(@(lambda)arrayfun(@(x)charfunBT(x),lambda),[-0.2,0.2]);
title('Characteristic function in BT for lambda in [-0.2,0.2]');
xlabel('\lambda');
ylabel('\Delta(\lambda)');
%% Plot Hopf bifurcation in two-parameter plane, including Takens-Bogdanov point
% and in $\beta-\delta-\omega$ space.
beta_hopf=getpar(hopf_branch_wbifs,ind.beta);
delta_hopf=getpar(hopf_branch_wbifs,ind.delta);
omega_hopf=[hopf_branch_wbifs.point.omega];
beta_bt=hbt.parameter(ind.beta);
delta_bt=hbt.parameter(ind.delta);
omega_bt=hbt.omega;
figure(3);clf;
subplot(1,2,1);
ax3l=gca;
phl=plot3(ax3l,beta_hopf,delta_hopf,omega_hopf,'b-',...
    beta_bt,delta_bt,omega_bt,'bo','linewidth',2);
grid on
view([-70,25]);
b2_lgtext_h={'Hopf','BT from Hopf'};
legend(ax3l,b2_lgtext_h,'location','northoutside');
xlabel(ax3l,'\beta');
ylabel(ax3l,'\delta');
zlabel(ax3l,'\omega');
subplot(1,2,2);
ax3r=gca;
ph=plot(ax3r,beta_hopf,delta_hopf,'b-',...
    beta_bt,delta_bt,'bo','linewidth',2);
grid on
xlabel(ax3r,'\beta');
ylabel(ax3r,'\delta');
%% Test functions for codimension-2 bifurcations along Hopf curve
% The test functions for 4 codimension-2 bifurcations are in output
% |hopftestfuncs| of the utility function HopfCodimension2:
%
% * L1 for generalized hopf (genh),
% * sign(omega) for Takens-Bogdanov (bt)
% * $\det(\Delta(0))$ for zero-Hopf interaction (zeho, touches zero for BT)
% * $\Re\lambda$ for complex eigenvalue closest to imaginary axis for
% Hopf-Hopf interaction (hoho).
%
figure(5);clf
plot(delta_hopf,[tanh(hopftestfuncs.genh(1,:));...
    hopftestfuncs.bt;...
    hopftestfuncs.zeho;...
    hopftestfuncs.hoho]);
grid on
legend({'tanh(genh)','BT','zero-Hopf','Hopf-Hopf'});
set(gca,'ylim',[-0.5,1.5]);
xlabel('\delta');
ylabel('\phi^h');
title('test functions along Hopf curve (Hopf-Hopf invisible)')
%% Fold bifurcation
% We do a standard continuation, using the starting index obtained from the
% flagged point indices.
fprintf('----- fold branch -----\n');
[fold_branch0,suc]=SetupFold(funcs,stst_branch_wbifs,br_getflags(stst_branch_wbifs,'fold'),...
    'contpar', [ind.delta,ind.beta], 'dir', ind.delta, 'step', 0.005, parameter_bd{:});
disp(['suc=',num2str(suc)]);
figure(2);
title(ax2,'Fold in beta-delta plane');
fold_branch0=br_contn(funcs,fold_branch0,300);
fold_branch0 = br_rvers(fold_branch0);
fold_branch0=br_contn(funcs,fold_branch0,300);
%% Detect special points along fold curve
fprintf('----- Codimension-two detection along fold branch -----\n');
[fold_branch_wbifs,foldtestfuncs,ind_foldcodim2,foldcodim2_type]=LocateSpecialPoints(funcs,fold_branch0);
nunst_fold=GetStability(fold_branch_wbifs,'exclude_trivial',true);
%% Insert fold in figure 3
beta_fold=getpar(fold_branch_wbifs,ind.beta);
delta_fold=getpar(fold_branch_wbifs,ind.delta);
omega_fold=0*beta_fold;
beta_bt2=bgetpar(fold_branch_wbifs,ind.beta,'BT');
delta_bt2=bgetpar(fold_branch_wbifs,ind.delta,'BT');
hold(ax3l,'on');
pf=plot3(ax3l,beta_fold,delta_fold,omega_fold,'k-',...
    beta_bt2,delta_bt2,0,'ks','linewidth',2);
b2_lgtext_f={'fold','BT from fold'};
legend(ax3l,[b2_lgtext_h,b2_lgtext_f],'location','northoutside');
hold(ax3r,'on');
plot(ax3r,beta_fold,delta_fold,'k-',...
    beta_bt2,delta_bt2,'ks','linewidth',2);
%% Test functions for codimension-2 bifurcations along fold curve
% The test functions for 4 codimension-2 bifurcations are in output
% |foldtestfuncs| of the utility function FoldCodimension2:
%
% * a for cusp (genh),
% * $[\det(\Delta(\lambda))]'$ in $\lambda=0$ for Takens-Bogdanov (bt)
% * $\Re\lambda$ for complex eigenvalue closest to imaginary axis for
% zero-Hopf interaction (zeho).
%
figure(5);clf
plot(delta_fold,[real(foldtestfuncs.cusp(1,:))/100;foldtestfuncs.bt;foldtestfuncs.zeho],'.');
grid on
legend({'cusp/100','BT','zero-Hopf'},'location','northwest');
set(gca,'ylim',[-1,1]);
xlabel('\delta');
ylabel('\phi^f');
title('test functions along fold curve (zero-Hopf invisible)')
%% Branch off to periodic orbit at some Hopf point, continue to large period
% We know that the orbit must be unstable close to the Hopf bifurcation
% from the Takens-Bogdanov normal form and the Lyapunov coefficients of the
% Hopf bifurcation.
psol_branch=SetupPsol(funcs,hopf_branch_wbifs,2,...
    'contpar',ind.delta,'degree',3,'intervals',50,parameter_bd{1:4},'max_step',[0,inf]);
[xm,ym]=df_measr(0,psol_branch);
ym.field='period';
ym.col=1;
ym.row=1;
psol_branch.method.continuation.plot_measure.x=xm;
psol_branch.method.continuation.plot_measure.y=ym;
figure(4);clf;
ax4=gca;
psol_branch=br_contn(funcs,psol_branch,25);
xlabel(ax4,'$\delta$','interpreter','latex');
ylabel(ax4,'period');
%% Convert point close to end of the psol_branch to homoclinic orbit
% and correct homoclinic orbit.Then refine and repeat correction.
hcli=p_tohcli(funcs,psol_branch.point(end-5));
figure(5);clf;
p_pplot(hcli);
xlabel('time/period');ylabel('x1,x2');
mhcli=df_mthod(funcs,'hcli');
mhcli.point.print_residual_info=1;
[hcli,suc]=p_correc(funcs,hcli,ind.delta,[],mhcli.point); % correct
disp(suc);
hcli2=p_remesh(hcli,3,50); % remesh it and
[hcli3,suc]=p_correc(funcs,hcli2,ind.beta,[],mhcli.point); % correct it again
disp(suc);
%% Continue branch of homoclinic orbits in two parameters
% Branches have to be created manually as described in hom_demo.
hcli_br=df_brnch(funcs,[ind.delta, ind.beta],'hcli');
hcli_br.point=hcli3;
hcli4=hcli3;
hcli4.parameter(ind.beta)=hcli4.parameter(ind.beta)-1e-4;
[hcli5,suc]=p_correc(funcs,hcli4,ind.delta,[],mhcli.point);
hcli_br.point(2)=hcli5;
hcli_br.parameter.max_bound=fold_branch_wbifs.parameter.max_bound;
hcli_br.parameter.min_bound=fold_branch_wbifs.parameter.min_bound;
hcli_br.parameter.max_step=[ind.beta,5e-3;ind.delta,5e-3];
hcli_br.method.point.print_residual_info=1;
figure(2);
hcli_br=br_contn(funcs,hcli_br,28);
hcli_br=br_rvers(hcli_br);
hcli_br=br_contn(funcs,hcli_br,40);
%% Add homoclinic orbit to two-parameter bifurcation diagram
beta_hcl=getpar(hcli_br,ind.beta);
delta_hcl=getpar(hcli_br,ind.delta);
phoml=plot3(ax3l,beta_hcl,delta_hcl,0*delta_hcl,'r-','linewidth',2);
phomr=plot(ax3r,beta_hcl,delta_hcl,'r-','linewidth',2);
legend(ax3l,[b2_lgtext_h,b2_lgtext_f,{'homoclinic'}],'location','northoutside');
%% Plot 2d diagram summary
figure(6);clf;ax6=gca;
hold(ax6,'on');
lg=Plot2dBranch(fold_branch_wbifs,'ax',ax6);
lg=Plot2dBranch(hopf_branch_wbifs,'ax',ax6,'oldlegend',lg);
lg=Plot2dBranch(hcli_br,'ax',ax6,'stability',0,'oldlegend',lg);
legend(lg{1},lg{2},'location','southeast');
xlabel('$\delta$','interpreter','latex');
ylabel('$\beta$','interpreter','latex');
