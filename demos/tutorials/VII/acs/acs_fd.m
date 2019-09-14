%%  Active control system
% 
% Demo illustrating how to branch off a Hopf Hopf bifurcation point
%
%% Differential equations for active control system with delayed feedback
%
% From:
% Yuting Ding, Jun Cao, Weihua Jiang. Double Hopf bifurcation in active 
% control system with delayed feedback: application to glue dosing 
% processes for particleboard. In Nonlinear Dynamics. 83, 1567--1576 (2016)
% 
% $$x'=(x+m)(1-x-m)-\frac{xy}{ay+x}$$ 
%
% $$y'=\delta y\left(\beta-\frac{y(t-\tau)}{x(t-\tau)}\right)$$
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id$
%
%% Clean workspace and add DDE-BifTool scripts to the MATLAB search path
clear;      % clear variables
close all;	% close figures
ddebiftoolpath='../../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
        strcat(ddebiftoolpath,'ddebiftool_extra_psol'),...
        strcat(ddebiftoolpath,'ddebiftool_extra_nmfm'),...
        strcat(ddebiftoolpath,'ddebiftool_utilities'));
%% Set parameter names
% This could also be loaded from mat file |'acs_parnames.mat'|.
parnames={'zeta','tau','tau_scaled'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Define the system
gu=0.1; gv=0.52; beta=0.1; % fixed parameters
acs_sys = @(xx,par) [...
    par(1,ind.tau,:).*xx(2,1,:);
    par(1,ind.tau,:).*(-xx(1,1,:)-gu*xx(1,2,:)-2*par(1,ind.zeta,:).*xx(2,1,:)-gv*xx(2,2,:)+beta*xx(1,2,:).^3)];
%% Set funcs structure
funcs=set_funcs('sys_rhs',acs_sys,'sys_tau',@()ind.tau_scaled,'x_vectorized',true,'p_vectorized',true);
%% Hopf-Hopf point
% Construct steady-state point
stst.kind='stst';
stst.x=[0;0];
stst.parameter(ind.zeta)=-0.016225;
stst.parameter(ind.tau)=5.89802;
stst.parameter(ind.tau_scaled)=1;
% Calculate stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1
%% Calculate coefficients of the parameter dependent normal form
hopf=p_tohopf(funcs,stst);
method=df_mthod(funcs,'hopf');
hopf.stability=p_stabil(funcs,hopf,method.stability);
hoho=p_tohoho(hopf);
unfolding_parameters=[ind.zeta, ind.tau];
hoho=nmfm_hoho(funcs,hoho,'free_pars',unfolding_parameters);
hoho.nmfm
%% Set bifurcation parameter range and step size bounds
brpars={'max_bound',[ind.tau 16],...
        'min_bound',[ind.tau 5],...
        'max_step', [ind.zeta 0.04; ind.tau 0.04]};
%% Continue Neimark-Sacker curves emanating from Hopf-Hopf point
[trfuncs,nsbr,suc]=C1branch_from_C2point(funcs,hoho,unfolding_parameters,...
    'codim2','hoho','codim1','TorusBifurcation',...
    brpars{:},'step',1e-01,'plot',0);
assert(all(suc(:)>0))
ntrsteps=186; [nsbr(1),suc]=br_contn(trfuncs,nsbr(1),ntrsteps);
assert(suc>0)
ntrsteps=61;  [nsbr(2),suc]=br_contn(trfuncs,nsbr(2),ntrsteps);
assert(suc>0)
%% Continue Hopf curve emanating from Hopf-Hopf point
[~,hbr,suc]=C1branch_from_C2point(funcs,hoho,unfolding_parameters,...
    'codim2','hoho','codim1','hopf',brpars{:},'step',1e-03,'plot',0);
assert(all(suc(:)>0))
nop=1000; [hbr(1),suc]=br_contn(funcs,hbr(1),nop);assert(suc>0)
hbr(1)=br_rvers(hbr(1));
[hbr(1),suc]=br_contn(funcs,hbr(1),nop); assert(suc>0)
nop=10; [hbr(2),suc]=br_contn(funcs,hbr(2),nop); assert(suc>0)
hbr(2)=br_rvers(hbr(2));
[hbr(2),suc]=br_contn(funcs,hbr(2),nop); assert(suc>0)
%% Plot the bifurcation curves using Plot2dBranch
% figure(1); clf; hold on;
% cm=colormap('lines');
% lg=Plot2dBranch(hbr(1), 'color',cm(1,:),'linewidth',1,'linestyle','-',...
%     'name','Hopf branches','stability',0);
% lg=Plot2dBranch(nsbr(1),'color',cm(2,:),'linewidth',1,'linestyle','-',...
%     'name','Neimark-Sacker branches','oldlegend',lg,'stability',0);
% lg=Plot2dBranch(nsbr(2),'color',cm(2,:),'linewidth',1,'linestyle','-',...
%     'name','Neimark-Sacker branches','oldlegend',lg,'stability',0);
% xlabel('\zeta');
% ylabel('\tau');
% title('Hopf curve and Neimark-Sacker curves emanating from the Hopf Hopf point')
%% Plot the bifurcation curves
figure(1); clf; hold on;
cm=colormap('lines');
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
hbr_pm = [getpars(hbr(1),ind.zeta); getpars(hbr(1),ind.tau)];
nsbr1_pm = [getpars(nsbr(1),ind.zeta); getpars(nsbr(1),ind.tau)];
nsbr2_pm = [getpars(nsbr(2),ind.zeta); getpars(nsbr(2),ind.tau)];
plot(hbr_pm(1,:),hbr_pm(2,:),'Color',cm(1,:),'DisplayName','Hopf branches')
plot(nsbr1_pm(1,:),nsbr1_pm(2,:),'Color',cm(2,:),'DisplayName','Neimark-Sacker branches')
plot(nsbr2_pm(1,:),nsbr2_pm(2,:),'Color',cm(2,:),'HandleVisibility','off')
plot(hoho.parameter(ind.zeta),hoho.parameter(ind.tau),'k.',...
    'MarkerSize',12,'DisplayName','Hopf Hopf point')
title('Hopf curve and Neimark-Sacker curves emanating from the Hopf Hopf point')
xlabel('$\zeta$','Interpreter','LaTex');
ylabel('$\tau$','Interpreter','LaTex');
legend; box on
%% Predictors for Neimark-Sacker and Hopf curves
[~,nsbr_pred]=C1branch_from_C2point(funcs,hoho,unfolding_parameters,...
    'codim2','hoho','codim1','TorusBifurcation',...
    'step',linspace(1e-03,2,40),'predictor',1);
[~,hbr_pred]=C1branch_from_C2point(funcs,hoho,unfolding_parameters,...
    'codim2','hoho','codim1','hopf',...
    'step',linspace(-2e-01,2e-01,30),'predictor',1);
%% Close-up near Hopf Hopf point in parameter space with predictors
figure(2); clf; hold on;
hbr1_pred_pm  = [getpars(hbr_pred(1), ind.zeta); getpars(hbr_pred(1), ind.tau)];
hbr2_pred_pm  = [getpars(hbr_pred(2), ind.zeta); getpars(hbr_pred(2), ind.tau)];
nsbr1_pred_pm = [getpars(nsbr_pred(1),ind.zeta); getpars(nsbr_pred(1),ind.tau)];
nsbr2_pred_pm = [getpars(nsbr_pred(2),ind.zeta); getpars(nsbr_pred(2),ind.tau)];
plot(hbr_pm(1,:),hbr_pm(2,:),'Color',cm(1,:),'DisplayName','Hopf branches')
plot(nsbr1_pm(1,:),nsbr1_pm(2,:),'Color',cm(2,:),'DisplayName','Neimark-Sacker branches')
plot(nsbr2_pm(1,:),nsbr2_pm(2,:),'Color',cm(2,:),'HandleVisibility','off')
plot(hbr1_pred_pm(1,:), hbr1_pred_pm(2,:) ,'.','Color',cm(1,:),'DisplayName','Hopf predictors')
plot(hbr2_pred_pm(1,:), hbr2_pred_pm(2,:) ,'.','Color',cm(1,:),'HandleVisibility','off')
plot(nsbr1_pred_pm(1,:),nsbr1_pred_pm(2,:),'.','Color',cm(2,:),'DisplayName','Neimark-Sacker predictors')
plot(nsbr2_pred_pm(1,:),nsbr2_pred_pm(2,:),'.','Color',cm(2,:),'HandleVisibility','off')
plot(hoho.parameter(ind.zeta),hoho.parameter(ind.tau),'k.',...
    'MarkerSize',12,'DisplayName','Hopf Hopf point')
title(['Close-up near Hopf Hopf point'...
      ' in parameter space with predictors'])
xlabel('$\zeta$','Interpreter','LaTex');
ylabel('$\tau$','Interpreter','LaTex');
axis([-0.0562    0.0109    5.7098    6.1757])
legend('Location','NorthWest'); box on
%%
% lg=Plot2dBranch(nsbr_pred(1),'color',cm(2,:),'name','Neimark-Sacker predictors',...
%     'Marker','.','linestyle','none','oldlegend',lg,'stability',0);
% lg=Plot2dBranch(nsbr_pred(2),'color',cm(2,:),'name','Neimark-Sacker predictors',...
%     'Marker','.','linestyle','none','oldlegend',lg,'stability',0);
% lg=Plot2dBranch(hbr_pred(1),'color',cm(1,:),'name','Hopf predictors',...
%     'Marker','.','linestyle','none','oldlegend',lg,'stability',0);
% lg=Plot2dBranch(hbr_pred(2),'color',cm(1,:),'name','Hopf predictors',...
%     'Marker','.','linestyle','none','oldlegend',lg,'stability',0);
% plot(hoho.parameter(ind.zeta),hoho.parameter(ind.tau),'k.',...
%     'MarkerSize',12,'DisplayName','Hopf Hopf point')
% xlabel('\zeta'); ylabel('\tau')
% legend('Location','NorthWest')
% axis([-0.0562    0.0109    5.7098    6.1757])
% box on
% title('Close-up near Hopf Hopf point in parameter space with predictors')
%% Compare computed and predicted periods
figure(3); clf; hold on;
% Plot computed periods on nsbr(1) and nsbr(2)
omegas1=arrayfun(@(p)p.period,nsbr(1).point);
omegas2=arrayfun(@(p)p.period,nsbr(2).point);
omegas1_pred=arrayfun(@(p)p.period,nsbr_pred(1).point);
omegas2_pred=arrayfun(@(p)p.period,nsbr_pred(2).point);
plot(getpars(nsbr(1),ind.zeta),omegas1,'Color',cm(1,:));
plot(getpars(nsbr(2),ind.zeta),omegas2,'Color',cm(2,:));
plot(getpars(nsbr_pred(1),ind.zeta),omegas1_pred,'.','Color',cm(1,:));
plot(getpars(nsbr_pred(2),ind.zeta),omegas2_pred,'.','Color',cm(2,:));
title('Compare computed and predicted periods')
xlabel('$\zeta$','Interpreter','LaTex');
ylabel('$\omega$','Interpreter','LaTex');
%% Detect codim-2 points on hopf_br1
[hopf_br_wbifs,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hbr(1));
% There are seven bifurcations detected on the Hopf branch hbr(1):
% four Hopf-Hopf bifurcations and three generalized Hopf bifurcations.
% However, by inspecting the parameters of the detected Hopf-Hopf points
% suggest that there are only two distinct Hopf Hopf points, 
% which are connected with the same Hopf branch. By comparing the location
% of the underlying points confirms the premise.
% We are therefore not interested in continuing the Hopf curves 
% emanating from these points. We are interested in continuing
% Neimark-Sacker curves emanating from the second Hopf-Hopf point
% and the limit point of cycle curves emanating from the three
% generalized Hopf points.
%% Bifurcation diagram using the function Plot2d
% figure(3);Plot2dBranch(hopf_br_wbifs);
% title('Bifurcation diagram')
% xlabel('$\zeta$','Interpreter','LaTex');
% ylabel('$\omega$','Interpreter','LaTex');
%% Subtract generalized Hopf points
genh_indices=hc2_indices(strcmp(hc2_types,'genh'));
genhpts=hopf_br_wbifs.point(genh_indices);
%% Initialize limit point of cycles emanating from the generalized Hopf points
[~,lpc_br1,suc]=C1branch_from_C2point(funcs,genhpts(1),...
    unfolding_parameters,'codim2','genh','codim1','POfold',...
    brpars{:},'step',1e-01,'plot',0);
assert(all(suc>0))
[poffuncs,lpc_br2,suc]=C1branch_from_C2point(funcs,genhpts(3),...
    unfolding_parameters,'codim2','genh','codim1','POfold',...
    brpars{:},'step',5e-03,'plot',0);
assert(all(suc>0))
%% Continue limit point of cycles emanating from the generalized Hopf points
ntrsteps=193; [lpc_br1,suc]=br_contn(poffuncs,lpc_br1,ntrsteps);
assert(suc>0)
ntrsteps=220; [lpc_br2,suc]=br_contn(poffuncs,lpc_br2,ntrsteps);
assert(suc>0)
%% Subtract second Hopf-Hopf point on the Hopf branch
hoho_indices=hc2_indices(strcmp(hc2_types,'hoho'));
hoho2=hopf_br_wbifs.point(hoho_indices(3));
%% Initialize and continue Neimark-Sacker curves emanating from Hopf-Hopf point
[trfuncs,nsbr34,suc]=C1branch_from_C2point(funcs,hoho2,unfolding_parameters,...
    'codim2','hoho','codim1','TorusBifurcation',...
    brpars{:},'step',5e-03,'plot',0);
assert(all(suc(:)>0))
% We only need to continue the second Neimark-Sacker curve since the first
% Neimark-Sacker curve emanating from second Hopf point is connected to the
% first Hopf point.
ntrsteps=200; [nsbr34(2),suc,fail,rjct]=br_contn(trfuncs,nsbr34(2),ntrsteps);
%% Bifurcation diagram in $(\zeta,\tau)$
figure(4); clf; hold on;
% Subtract paramater values from the branches
lpc_br1_pm= [getpars(lpc_br1,   ind.zeta); getpars(lpc_br1,  ind.tau)];
lpc_br2_pm= [getpars(lpc_br2,   ind.zeta); getpars(lpc_br2,  ind.tau)];
nsbr3_pm  = [getpars(nsbr34(2), ind.zeta); getpars(nsbr34(2),ind.tau)];
% Plot curves
plot(hbr_pm(1,:),    hbr_pm(2,:),    'Color',cm(1,:),'DisplayName','Hopf branches');
plot(nsbr1_pm(1,:),  nsbr1_pm(2,:),  'Color',cm(2,:),'DisplayName','Neimark-Sacker branches');
plot(nsbr2_pm(1,:),  nsbr2_pm(2,:),  'Color',cm(2,:),'HandleVisibility','off');
plot(nsbr3_pm(1,:),  nsbr3_pm(2,:),  'Color',cm(2,:),'HandleVisibility','off');
plot(lpc_br1_pm(1,:),lpc_br1_pm(2,:),'Color',cm(3,:),'DisplayName','Generalized Hopf');
plot(lpc_br2_pm(1,:),lpc_br2_pm(2,:),'Color',cm(3,:),'HandleVisibility','off');
% Add bifurcation points
plot(hoho.parameter(ind.zeta), hoho.parameter(ind.tau),...
    'k.','MarkerSize',8,'DisplayName','Hopf Hopf point')
plot(hoho2.parameter(ind.zeta),hoho2.parameter(ind.tau),...
    'k.','MarkerSize',8,'HandleVisibility','off')
plot(genhpts(1).parameter(ind.zeta),genhpts(1).parameter(ind.tau),...
    'b.','MarkerSize',8,'DisplayName','Generelized Hopf point')
plot(genhpts(2).parameter(ind.zeta),genhpts(2).parameter(ind.tau),...
    'b.','MarkerSize',8,'HandleVisibility','off')
plot(genhpts(3).parameter(ind.zeta),genhpts(3).parameter(ind.tau),...
    'b.','MarkerSize',8,'HandleVisibility','off')
title('Bifurcation diagram in (\zeta,\tau)')
xlabel('$\zeta$','Interpreter','LaTex')
ylabel('$\tau$','Interpreter','LaTex')
axis([-0.3    0.4    5    16])
legend()
%% Plot paramater points which are used for simulation near the first Hopf Hopf point
figure(5); clf; hold('on')
% Connect ns1_br and ns2_br to hoho for nice graphic
nsbr(1).point=[nsbr(1).point(1) nsbr(1).point];
nsbr(1).point(1).parameter(1:3)=hoho.parameter;
nsbr(2).point=[nsbr(2).point(1) nsbr(2).point];
nsbr(2).point(1).parameter(1:3)=hoho.parameter;
% Subtract paramater values from the branches
hbr2_pm=[getpars(hbr(2),ind.zeta); getpars(hbr(2),ind.tau)];
nsbr1_pm=[getpars(nsbr(1),ind.zeta); getpars(nsbr(1),ind.tau)];
nsbr2_pm=[getpars(nsbr(2),ind.zeta); getpars(nsbr(2),ind.tau)];
% Plot bifurcation curves
indices=find(hbr_pm(1,:)>-0.0170 & hbr_pm(1,:)<-0.015);
plot(hbr_pm(1,indices),hbr_pm(2,indices),'Color',cm(1,:),'HandleVisibility','off')
plot(hbr2_pm(1,:),hbr2_pm(2,:),'Color',cm(1,:),'HandleVisibility','off')
plot(nsbr1_pm(1,:),nsbr1_pm(2,:),'Color',cm(2,:),'HandleVisibility','off')
plot(nsbr2_pm(1,:),nsbr2_pm(2,:),'Color',cm(2,:),'HandleVisibility','off')
% Add Hopf Hopf point
plot(hoho.parameter(ind.zeta),hoho.parameter(ind.tau),...
    'k.','MarkerSize',8,'DisplayName','Hopf Hopf point')
% Calculate parameters used for simulation
tau=5.901783308978358;
zeta=-0.015485728828307;
zeta1 = zeta;
zeta2 = zeta-0.0002;
zeta3 = zeta-0.0004;
zeta4 = zeta-0.000445;
zeta5 = zeta-0.0004496;
% Add to plot
plot(zeta1,tau,'.','DisplayName','$\zeta_1$')
plot(zeta2,tau,'.','DisplayName','$\zeta_2$')
plot(zeta3,tau,'.','DisplayName','$\zeta_3$')
plot(zeta4,tau,'.','DisplayName','$\zeta_4$')
plot(zeta5,tau,'.','DisplayName','$\zeta_5$')
% Add legend, axis labels, plot title and set range
title('Simulation points near Hopf Hopf point')
xlabel('$\zeta$','Interpreter','LaTex','FontSize',9)
ylabel('$\tau$','Interpreter','LaTex','FontSize',9)
axis([-0.0170   -0.0152    5.8968    5.9023])
legend('Location','NorthWest','Interpreter','LaTeX')
box on