%%  Rose-Hindmarsh model with time delay
% 
% Demo illustrating how to branch off a fold-Hopf bifurcation point
%
%% Differential equations for Rose-Hindmarsh model with time delay
%
% From:
% S. Ma and Z. Feng,
% Fold-Hopf bifurcation of the Rose-Hindmarsh model with time delay,
% Appl. Math. Comput., 219 (2013), pp. 10073--10081
% https://www.sciencedirect.com/science/article/abs/pii/S0096300313004141
% 
% $$\dot{x}=y - ax^3 + bx^2(t-tau) - cz + I_{app},$$
%
% $$\dot{y}=c - dx^2-y,$$
%
% $$\dot{z]=r(S(x-\chi)-z).$$
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
parnames={'Iapp','S','r','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Define the system
a=1.0; b=3.0; c=1.0; d=5.0; chi=-1.6;
RH_sys = @(xx,par) [...
    xx(2,1,:)-a*xx(1,1,:).^3+b*xx(1,2,:).^2-c*xx(3,1,:)+par(1,ind.Iapp,:);
    c-d*xx(1,1,:).^2-xx(2,1,:);
    par(1,ind.r,:).*(par(1,ind.S,:).*(xx(1,1,:)-chi)-xx(3,1,:))];
%% Set funcs structure
funcs=set_funcs('sys_rhs',RH_sys,'sys_tau',@()ind.tau,...
    'x_vectorized',true,'p_vectorized',true);
%% Construct fold-Hopf point
a=1.0; b=3.0; c=1.0; d=5.0; chi=-1.6;
r=1.0e-03;
S=-0.57452592;
[xstar,Ip,tau]=bifurcationvalues(a,b,c,d,chi,r,S);
% % Construct steady-state point
stst=dde_stst_create('x',[xstar; c-d*xstar^2; S*(xstar-chi)]);
stst.parameter([ind.Iapp ind.S ind.r ind.tau])=[Ip S r tau];
% Calculate stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1(1:5)
% Calculate normal coefficients
hopf=p_tohopf(funcs,stst);
zeho=p_tozeho(hopf);
unfolding_pars=[ind.Iapp, ind.S];
zeho=nmfm_zeho(funcs,zeho,'free_pars',unfolding_pars);
zeho.nmfm
%% Set bifurcation parameter range and step size bounds
brpars={'min_bound',[ind.Iapp -20; ind.S -12],...
        'max_bound',[ind.Iapp  10; ind.S  5],...
        'max_step', [ind.Iapp 4e-02; ind.S 4e-02]};
%% Continue Neimark-Sacker curve emanating from fold-Hopf point
[trfuncs,nsbr,~]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2',zeho.kind,'codim1','TorusBifurcation',...
    brpars{:},'step',1e-03,'plot',0);
ntrsteps=1000; [nsbr,suc]=br_contn(trfuncs,nsbr,ntrsteps);
% Calculate stability on the Neimark-Sacker branch
% [nunst,dom,triv_defect,nsbr.point]=GetStability(nsbr,'funcs',trfuncs);
% nunst'
%% Continue Hopf curve emanating from fold-Hopf point
[~,hbr,suc]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2',zeho.kind,'codim1','hopf',brpars{:},'step',1e-03,'plot',0);
assert(all(suc(:)>0))
nop=1000; [hbr,suc]=br_contn(funcs,hbr,nop); assert(suc>0)
hbr=br_rvers(hbr);
[hbr,suc]=br_contn(funcs,hbr,nop); assert(suc>0)
%% Continue fold curve emanating from fold-Hopf point
[~,fbr,suc]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2',zeho.kind,'codim1','fold',brpars{:},...
    'step',1e-03,'plot',0);
nop=1000; [fbr,suc]=br_contn(funcs,fbr,nop);
fbr=br_rvers(fbr);
[fbr,suc]=br_contn(funcs,fbr,nop);
%% Detect special points on the Hopf branch
[hbr_wbifs,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hbr);
Ip=arrayfun(@(x)x.parameter(ind.Iapp),hbr_wbifs.point);
figure(6); clf; hold on
plot(Ip,hopftests.zeho(1,:),'.-')
plot(Ip,hopftests.genh(1,:),'.-');
plot(Ip,zeros(size(Ip)));
xlabel('$I_{app}$','Interpreter','LaTex');
ylabel('First Lyapunov coefficient (L1)')
title('Criticality along Hopf bifurcation curve')
legend('zeho-Hopf','generalized Hopf','zero')
%% Plot the bifurcation curves
figure(1); clf; hold on;
cm=colormap('lines');
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
hbr_pm    = [getpars(hbr,ind.Iapp); getpars(hbr,ind.S)];
nsbr1_pm  = [getpars(nsbr,ind.Iapp); getpars(nsbr,ind.S)];
fbr_pm    = [getpars(fbr,ind.Iapp);  getpars(fbr,ind.S)];
plot(hbr_pm(1,:),hbr_pm(2,:),'Color',cm(1,:),'DisplayName','Hopf branch');
plot(nsbr1_pm(1,:),nsbr1_pm(2,:),'Color',cm(3,:),...
    'DisplayName','Neimark-Sacker branch');
plot(fbr_pm(1,:),fbr_pm(2,:),'Color',cm(5,:),'DisplayName','Fold branch');
plot(zeho.parameter(ind.Iapp),zeho.parameter(ind.S),'k.',...
    'MarkerSize',12,'DisplayName','fold-Hopf point')
title('Bifurcation curves emanating from the fold-Hopf point')
axis([-20 10 -12 5])
xlabel('$I_{app}$','Interpreter','LaTex');
ylabel('$S$','Interpreter','LaTex');
legend('Location','NorthWest'); box on
%% Plot the bifurcation curves (with hopf criticality)
figure(1); clf; hold on;
L1s=hopftests.genh(1,:);
hbrsub.point=hbr_wbifs.point(L1s>0);
hbrsuper.point=hbr_wbifs.point(L1s<0);
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
hbrsub_pm = [getpars(hbrsub,ind.Iapp); getpars(hbrsub,ind.S)];
hbrsup_pm = [getpars(hbrsuper,ind.Iapp); getpars(hbrsuper,ind.S)];
nsbr1_pm  = [getpars(nsbr,ind.Iapp); getpars(nsbr,ind.S)];
fbr_pm    = [getpars(fbr,ind.Iapp);  getpars(fbr,ind.S)];
plot(hbrsub_pm(1,:),hbrsub_pm(2,:),'Color',cm(1,:),...
    'DisplayName','subcritical Hopf branch');
plot(hbrsup_pm(1,:),hbrsup_pm(2,:),'Color',cm(2,:),...
    'DisplayName','supercritical Hopf branch');
plot(nsbr1_pm(1,:),nsbr1_pm(2,:),'Color',cm(3,:),...
    'DisplayName','Neimark-Sacker branch');
plot(fbr_pm(1,:),fbr_pm(2,:),'Color',cm(5,:),'DisplayName','Fold branch');
plot(zeho.parameter(ind.Iapp),zeho.parameter(ind.S),'k.',...
    'MarkerSize',12,'DisplayName','fold-Hopf point')
title('Bifurcation curves emanating from the fold-Hopf point')
axis([-20 10 -12 5])
xlabel('$I_{app}$','Interpreter','LaTex');
ylabel('$S$','Interpreter','LaTex');
legend('Location','NorthWest'); box on
%% Predictors for Neimark-Sacker and Hopf curves
[~,nsbr_pred]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2','zeho','codim1','TorusBifurcation',...
    'step',linspace(1e-03,2e-01,100),'predictor',true);
[~,hbr_pred]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2','zeho','codim1','hopf',brpars{:},...
    'step',linspace(1e-03,1e-02),'predictor',1);
%% Plot comparing computed and predicted Neimark-Sacker curve
figure(2); clf; hold on;
nsbr2_pm_pred = [getpars(nsbr_pred,ind.Iapp); getpars(nsbr_pred,ind.S)];
plot(nsbr1_pm(1,:),nsbr1_pm(2,:),'Color',cm(2,:),...
    'DisplayName','Neimark-Sacker branch');
plot(nsbr2_pm_pred(1,:),nsbr2_pm_pred(2,:),'.','Color',cm(2,:),...
    'DisplayName','Neimark-Sacker predictor');
title('Neimark-Sacker curve emanating from the fold-Hopf point')
axis([-20 2 -12 2])
plot(zeho.parameter(ind.Iapp),zeho.parameter(ind.S),'k.',...
    'MarkerSize',12,'DisplayName','fold-Hopf point')
xlabel('$I_{app}$','Interpreter','LaTex');
ylabel('$S$','Interpreter','LaTex');
legend('Location','NorthWest')
box on
%% Plot comparing computed and predicted periodic orbits
figure(3); clf; hold on;
for i=1:11
   plot3(nsbr.point(i).profile(1,:),nsbr.point(i).profile(2,:),...
       nsbr.point(i).profile(3,:),'Color',cm(1,:));
end
for i=1:8
   plot3(nsbr_pred.point(i).profile(1,:),nsbr_pred.point(i).profile(2,:),...
       nsbr_pred.point(i).profile(3,:),'Color',cm(2,:));
end
title('Comparison between computed and predicted periodic orbits')
view(3)
xlabel('x'); ylabel('y'); zlabel('z')
%% Detect special points on the Hopf branch
[hbr_wbifs,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hbr);
Ip=arrayfun(@(x)x.parameter(ind.Iapp),hbr_wbifs.point);
figure(4); clf;
plot(Ip,hopftests.genh(1,:),'.-',Ip,zeros(size(Ip)));
xlabel('$I_{app}$','Interpreter','LaTex');
ylabel('First Lyapunov coefficient (L1)')
title('Criticality along Hopf bifurcation curve')
axis([-20 7 -15 15]);
%% Bifurcation diagram in (x^{*},S) space as in the article
figure(5); clf; hold on
plot(arrayfun(@(p) p.x(1), hbr.point), getpars(hbr,ind.S))
plot(arrayfun(@(p) p.x(1), fbr.point), getpars(fbr,ind.S))
axis([0.06, 0.22, -1,-0.2])
xlabel('x^*'); ylabel('S')
legend('Hopf branch','Fold branch')
%% Different parameters with stable torus
r=1.4; S=-8;
[xstar,Ip,tau]=bifurcationvalues(a,b,c,d,chi,r,S);
% Construct steady-state point
stst=dde_stst_create('x',[xstar; c-d*xstar^2; S*(xstar-chi)]);
stst.parameter([ind.Iapp ind.S ind.r ind.tau])=[Ip S r tau];
% Calculate stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1(1:5)
% Calculate normal coefficients
hopf=p_tohopf(funcs,stst);
zeho=p_tozeho(hopf);
unfolding_pars=[ind.Iapp, ind.S];
zeho=nmfm_zeho(funcs,zeho,'free_pars',unfolding_pars);
fprintf('s=%f, theta=%f, e=%f\n',...
	zeho.nmfm.s,zeho.nmfm.theta,zeho.nmfm.e)
%% Set bifurcation parameter range and step size bounds
brpars={'min_bound',[ind.Iapp -20; ind.S -9],...
        'max_bound',[ind.Iapp -18; ind.S -7],...
        'max_step', [ind.Iapp 0.04; ind.S 0.04]};
%% Continue Neimark-Sacker curve emanating from fold-Hopf point
[trfuncs,nsbr,suc]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2','zeho','codim1','TorusBifurcation',...
    brpars{:},'step',1e-03,'plot',0);
assert(all(suc(:)>0))
ntrsteps=1000; [nsbr,suc]=br_contn(trfuncs,nsbr,ntrsteps);
%% Continue Hopf curve emanating from fold-Hopf point
[~,hbr,suc]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2','zeho','codim1','hopf',brpars{:},'step',1e-05,'plot',0);
assert(all(suc(:)>0))
nop=1000; [hbr,suc]=br_contn(funcs,hbr,nop);
hbr=br_rvers(hbr);
[hbr,suc]=br_contn(funcs,hbr,nop); assert(suc>0)
%% Continue fold curve emanating from fold-Hopf point
[~,fbr,suc]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2','zeho','codim1','fold',brpars{:},'step',1e-03,'plot',0);
assert(all(suc(:)>0))
nop=100; [fbr,suc]=br_contn(funcs,fbr,nop); assert(suc>0)
fbr=br_rvers(fbr);
[fbr,suc]=br_contn(funcs,fbr,nop); assert(suc>0)
%% Detect special points on the Hopf branch
[hbr_wbifs,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hbr);
Ip=arrayfun(@(x)x.parameter(ind.Iapp),hbr_wbifs.point);
figure(6); clf; hold on
plot(Ip,hopftests.zeho(1,:),'.-')
plot(Ip,hopftests.genh(1,:),'.-');
plot(Ip,zeros(size(Ip)));
xlabel('$I_{app}$','Interpreter','LaTex');
ylabel('First Lyapunov coefficient (L1)')
title('Criticality along Hopf bifurcation curve')
legend('zeho-Hopf','generalized Hopf','zero')
%% Plot the bifurcation curves
figure(7); clf; hold on;
cm=colormap('lines');
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
L1s=hopftests.genh(1,:);
hbrsub.point=hbr.point(L1s>0);
hbrsup.point=hbr.point(L1s<0);
hbrsub_pm   = [getpars(hbrsub,ind.Iapp); getpars(hbrsub,ind.S)];
hbrsup_pm = [getpars(hbrsup,ind.Iapp); getpars(hbrsup,ind.S)];
nsbr1_pm = [getpars(nsbr,ind.Iapp); getpars(nsbr,ind.S)];
fbr_pm = [getpars(fbr,ind.Iapp); getpars(fbr,ind.S)];
plot(hbrsub_pm(1,:),hbrsub_pm(2,:),'Color',cm(1,:),...
    'DisplayName','subcritical Hopf branch');
plot(hbrsup_pm(1,:),hbrsup_pm(2,:),'Color',cm(2,:),...
    'DisplayName','supercritical Hopf branch');
plot(nsbr1_pm(1,:),nsbr1_pm(2,:),'Color',cm(3,:),'DisplayName','Neimark-Sacker branch');
% plot(fbr_pm(1,:),fbr_pm(2,:),'Color',cm(3,:),'DisplayName','Fold branch');
plot(zeho.parameter(ind.Iapp),zeho.parameter(ind.S),'k.',...
    'MarkerSize',12,'DisplayName','fold-Hopf point')
title('Bifurcation curves emanating from the fold-Hopf point')
axis([-19.0193  -18.7128   -8.0477   -7.9587])
xlabel('$I_{app}$','Interpreter','LaTex');
ylabel('$S$','Interpreter','LaTex');
legend('Location','NorthWest')
box on
%% Predictors for Neimark-Sacker and Hopf curves
[~,nsbr_pred]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2','zeho','codim1','TorusBifurcation',...
    'step',linspace(1e-03,2.2e-01,40),'predictor',1);
[~,hbrsub_pred]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2','zeho','codim1','hopf',...
    'step',linspace(0,1e-03,20),'predictor',1);
[~,hbrsup_pred]=C1branch_from_C2point(funcs,zeho,unfolding_pars,...
    'codim2','zeho','codim1','hopf',...
    'step',linspace(-1e-03,0,20),'predictor',1);
%% Plot comparing computed and predicted Neimark-Sacker curve
figure(8); clf; hold on;
nsbr2_pm_pred  = [getpars(nsbr_pred,ind.Iapp);   getpars(nsbr_pred,ind.S)];
hbrsub_pm_pred = [getpars(hbrsub_pred,ind.Iapp); getpars(hbrsub_pred,ind.S)];
hbrsup_pm_pred = [getpars(hbrsup_pred,ind.Iapp); getpars(hbrsup_pred,ind.S)];
plot(hbrsub_pm(1,:),hbrsub_pm(2,:),'Color',cm(1,:),...
    'DisplayName','subcritical Hopf branch');
plot(hbrsup_pm(1,:),hbrsup_pm(2,:),'Color',cm(2,:),...
    'DisplayName','supercritical Hopf branch');
plot(hbrsub_pm_pred(1,:),hbrsub_pm_pred(2,:),'.','Color',cm(1,:),...
    'DisplayName','subcritical Hopf predictor');
plot(hbrsup_pm_pred(1,:),hbrsup_pm_pred(2,:),'.','Color',cm(2,:),...
    'DisplayName','supercritical Hopf predictor');
plot(nsbr1_pm(1,:),nsbr1_pm(2,:),'Color',cm(3,:),...
    'DisplayName','Neimark-Sacker branch');
plot(nsbr2_pm_pred(1,:),nsbr2_pm_pred(2,:),'.','Color',cm(3,:),...
    'DisplayName','Neimark-Sacker predictor');
plot(zeho.parameter(ind.Iapp),zeho.parameter(ind.S),'k.',...
    'MarkerSize',12,'DisplayName','fold-Hopf point')
title('Neimark-Sacker curve emanating from the fold-Hopf point')
axis([-19.0193  -18.7128   -8.0477   -7.9587])
xlabel('$I_{app}$','Interpreter','LaTex');
ylabel('$S$','Interpreter','LaTex');
text(-18.784,-7.981,  'I','FontSize',14); % stable period orbit
text(-18.905,-8.032, 'II','FontSize',14); % stable 2d torus
legend('Location','NorthWest'); box on
%% Plot comparing computed and predicted periodic orbits
figure(9); clf; hold on;
for i=1:14
   plot3(nsbr.point(i).profile(1,:),nsbr.point(i).profile(2,:),...
         nsbr.point(i).profile(3,:),'Color',cm(1,:));
end
for i=1:9
   plot3(nsbr_pred.point(i).profile(1,:),nsbr_pred.point(i).profile(2,:),...
         nsbr_pred.point(i).profile(3,:),'Color',cm(2,:));
end
title('Comparison between computed and predicted periodic orbits')
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','LaTex');
zlabel('$z$','Interpreter','LaTex')
view(3); grid on
%% Compare computed and predicted periods
figure(10); clf; hold on;
% Plot computed periods on nsbr(1) and nsbr(2)
omegas1=arrayfun(@(p)p.period,nsbr.point);
omegas1_pred=arrayfun(@(p)p.period,nsbr_pred.point);
plot(getpars(nsbr,ind.Iapp),omegas1,'Color',cm(1,:));
plot(getpars(nsbr_pred,ind.Iapp),omegas1_pred,'.','Color',cm(1,:));
title('Compare computed and predicted periods')
xlabel('$I_{app}$','Interpreter','LaTex')
ylabel('$\omega$','Interpreter','LaTex')
axis([-19.1016  -18.6500    1.1144    1.1220])
legend('Computed period','Predicted period','Location','SouthEast'); box on