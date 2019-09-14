%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id$
%% clean workspace and add DDE-BifTool scripts to search path
clear;      % clear variables
close all;	% close figures
ddebiftoolpath='../../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_psol'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_nmfm'),...
    strcat(ddebiftoolpath,'ddebiftool_utilities'));
%% Set parameter names
% This could also be loaded from mat file |'vdpo_parnames.mat'|.
parnames={'mu1','mu2','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Set funcs structure
% We load the precalculated multilinear forms. These have been generated
% with the file gen_sym_vdpo.m.
funcs=set_symfuncs(@sym_vdpo_mf,'sys_tau',@()ind.tau);
%% Construct transcritical-Hopf bifucation point
stst=dde_stst_create('x',[0;0],'parameter',[0 0 1]);
% Calculate stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1
%% Coefficients of the parameter dependent normal form
ht=p_tohopf(funcs,stst);
% ht=dde_hopf_clear(ht);
ht=p_tozeho(ht);
unfolding_pars=[ind.mu1, ind.mu2];
ht=nmfm_zeho(funcs,ht,'transcritical',1,'free_pars',unfolding_pars);
ht.nmfm
%% Set bifurcation parameter range and step size bounds
brpars={'max_bound',[ind.mu1 0.0139; ind.mu2 0.0094],...
    'min_bound',[ind.mu1 -0.0190; ind.mu2 -0.0069 ],...
    'max_step' ,[ind.mu1  1.0e-02; ind.mu2 1.0e-02]};
%% Continue Neimark-Sacker curves emanating from transcritical-Hopf point
[trfuncs,nsbr,suc]=C1branch_from_C2point(funcs,ht,unfolding_pars,...
    'codim2','zeho','codim1','TorusBifurcation',...
    'step',1e-04,'plot',0,brpars{:});
assert(all(suc(:)>0))
ntrsteps=27; [nsbr(1),suc]=br_contn(trfuncs,nsbr(1),ntrsteps); assert(suc>0)
ntrsteps=30; [nsbr(2),suc]=br_contn(trfuncs,nsbr(2),ntrsteps); assert(suc>0)
%% Continue Hopf curves emanating from fold-Hopf point
[~,hbr,suc]=C1branch_from_C2point(funcs,ht,unfolding_pars,...
    'codim2','zeho','codim1','hopf',brpars{:},'step',1e-05,'plot',0);
assert(all(suc(:)>0))
for i=2:-1:1
    nop=1000; hbr(i)=br_contn(funcs,hbr(i),nop);
    hbr(i)=br_rvers(hbr(i));
    hbr(i)=br_contn(funcs,hbr(i),nop);
end
%% Continue transcritical curve emanating from fold-Hopf point
[~,tcbr]=C1branch_from_C2point(funcs,ht,unfolding_pars,...
    'codim2','zeho','codim1','fold',brpars{1:4},....
    'step',linspace(-8.0e-03,8.0e-03,10),'plot',0);
%% Detect special points on the Hopf branches
for i=2:-1:1
    [hbr_wbifs(i),hopftests(i),hc2_indices,hc2_types]=...
        LocateSpecialPoints(funcs,hbr(i));
    al{i}=arrayfun(@(x)x.parameter(ind.mu1),hbr_wbifs(i).point);
    figure(i); clf;
    plot(al{i},hopftests(i).zeho(1,:),'.-',al{i},zeros(size(al{i})));
    xlabel('$\mu_1$','Interpreter','LaTex');
    ylabel('First Lyapunov coefficient (L1)')
    title('Criticality along Hopf bifurcation curve')
end
%% Plot the bifurcation curves
figure(3); clf; hold on;
cm=colormap('lines');
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
plot(ht.parameter(ind.mu1),ht.parameter(ind.mu2),'k.','MarkerSize',12)
tcbr_pm = [getpars(tcbr,ind.mu1); getpars(tcbr,ind.mu2)];
plot(tcbr_pm(1,:),tcbr_pm(2,:),'Color',cm(5,:));
for i=2:-1:1
    L1s=hopftests(i).zeho(1,:);
    hbrsub(i).point=hbr(i).point(L1s>0);
    hbrsup(i).point=hbr(i).point(L1s<0);
    hbrsub_pm{i} = [getpars(hbrsub(i),ind.mu1); getpars(hbrsub(i),ind.mu2)];
    hbrsup_pm{i} = [getpars(hbrsup(i),ind.mu1); getpars(hbrsup(i),ind.mu2)];
    nsbr_pm{i} = [getpars(nsbr(i),ind.mu1); getpars(nsbr(i),ind.mu2)];
    plot(nsbr_pm{i}(1,:),nsbr_pm{i}(2,:),'Color',cm(3,:));
    plot(hbrsup_pm{i}(1,:),hbrsup_pm{i}(2,:),'Color',cm(2,:));
    plot(hbrsub_pm{i}(1,:),hbrsub_pm{i}(2,:),'Color',cm(1,:));
end
legend({'transcritical-Hopf','transcritical curve','Neimark-Sacker branch',...
    'superscritical Hopf branch','subcritical Hopf branch'})
title('Bifurcation curves emanating from the fold-Hopf point')
axis([-0.0190    0.0139   -0.0069    0.0094])
xlabel('$\mu_1$','Interpreter','LaTex');
ylabel('$\mu_2$','Interpreter','LaTex');
box on
% Reverse the stacking order of the graphics
chi=get(gca,'Children'); set(gca,'Children',flipud(chi));
%% Predictors for Neimark-Sacker, Hopf and transcritical curves
[trfuncs,nsbr_pred]=C1branch_from_C2point(funcs,ht,unfolding_pars,...
    'codim2','zeho','codim1','TorusBifurcation',...
    'step',linspace(1e-05,4e-02,20),'predictor',true);
[~,tcbr_pred]=C1branch_from_C2point(funcs,ht,unfolding_pars,...
    'codim2','zeho','codim1','fold',...
    'step',linspace(-8.0e-03,8.0e-03,10),'predictor',true);
[~,hbrsub_pred]=C1branch_from_C2point(funcs,ht,unfolding_pars,...
    'codim2','zeho','codim1','hopf',brpars{:},...
    'step',linspace(-1e-01,0,80),'predictor',1);
[~,hbrsup_pred]=C1branch_from_C2point(funcs,ht,unfolding_pars,...
    'codim2','zeho','codim1','hopf',brpars{:},...
    'step',linspace(0,1e-01,40),'predictor',1);
% Correct super- and subcritical hopf branches for the first curve
tempbr=hbrsub_pred(1);
hbrsub_pred(1)=hbrsup_pred(1);
hbrsup_pred(1)=tempbr;
%% Plot comparing computed and predicted curves
figure(4); clf; hold on;
plot(ht.parameter(ind.mu1),ht.parameter(ind.mu2),'k.','MarkerSize',12)
tcbr_pm_pred = [getpars(tcbr_pred,ind.mu1); getpars(tcbr_pred,ind.mu2)];
plot(tcbr_pm_pred(1,:),tcbr_pm_pred(2,:),'.','Color',cm(5,:));
plot(tcbr_pm(1,:),tcbr_pm(2,:),'Color',cm(5,:));
for i=2:-1:1
    nsbr_pm_pred{i} = [getpars(nsbr_pred(i),ind.mu1); ...
        getpars(nsbr_pred(i),ind.mu2)];
    hbrsub_pm_pred{i}  = [getpars(hbrsub_pred(i),ind.mu1); ...
        getpars(hbrsub_pred(i),ind.mu2)];
    hbrsup_pm_pred{i} = [getpars(hbrsup_pred(i),ind.mu1); ...
        getpars(hbrsup_pred(i),ind.mu2)];
    plot(nsbr_pm_pred{i}(1,:),nsbr_pm_pred{i}(2,:),'.','Color',cm(3,:));
    plot(hbrsub_pm_pred{i}(1,:),hbrsub_pm_pred{i}(2,:),'.','Color',cm(2,:));
    plot(hbrsup_pm_pred{i}(1,:),hbrsup_pm_pred{i}(2,:),...
        '.','Color',cm(1,:));
    plot(nsbr_pm{i}(1,:),nsbr_pm{i}(2,:),'Color',cm(3,:));
    plot(hbrsup_pm{i}(1,:),hbrsup_pm{i}(2,:),'Color',cm(2,:));
    plot(hbrsub_pm{i}(1,:),hbrsub_pm{i}(2,:),'Color',cm(1,:));
end
legend({'transcritical Hopf','transcritical predictor',...
    'transcritical branch','Neimark-Sacker predictor',...
    'supercritical Hopf predictor','subcritical Hopf predictor',...
    'Neimark-Sacker branches','supercritical Hopf branches',...
    'subscritical Hopf branches'})
title('Neimark-Sacker curve emanating from the transcrical-Hopf point')
axis([-0.0190    0.0139   -0.0069    0.0094])
xlabel('$\mu_1$','Interpreter','LaTex')
ylabel('$\mu_2$','Interpreter','LaTex')
text(-0.0129,0.0038,'I');
text(-0.0093,0.007,'II');
text(0.008,-0.0052,'III');
text(0.0115,-0.003,'IV');
legend('Location','NorthEast'); box on
% reverse the stacking order of the graphics
chi=get(gca, 'Children');set(gca, 'Children',flipud(chi));
%% Plot comparing computed and predicted periodic orbits
figure(5); clf; hold on;
genpars=@(br,i)ones(1,length(br.point(1).profile(1,:))).*br.point(i).parameter(ind.mu1);
for i=1:23
    plot3(genpars(nsbr(1),i),nsbr(1).point(i).profile(1,:),...
        nsbr(1).point(i).profile(2,:),'Color',cm(1,:));
end
for i=1:12
    plot3(genpars(nsbr_pred(1),i),nsbr_pred(1).point(i).profile(1,:),...
        nsbr_pred(1).point(i).profile(2,:),'Color',cm(2,:));
end
for i=1:20
    plot3(genpars(nsbr(2),i),nsbr(2).point(i).profile(1,:),...
        nsbr(2).point(i).profile(2,:),'Color',cm(1,:));
end
for i=1:12
    plot3(genpars(nsbr_pred(2),i),nsbr_pred(2).point(i).profile(1,:),...
        nsbr_pred(2).point(i).profile(2,:),'Color',cm(2,:));
end
title('Comparison between computed and predicted periodic orbits')
xlabel('$\mu_1$','Interpreter','LaTex');
ylabel('$x$','Interpreter','LaTex');
zlabel('$y$','Interpreter','LaTex');
view(3)
%% Compare computed and predicted periods
figure(6); clf; hold on;
% plot computed periods on nsbr(1) and nsbr(2)
omegas1=arrayfun(@(p)p.period,nsbr(1).point);
omegas2=arrayfun(@(p)p.period,nsbr(2).point);
omegas1_pred=arrayfun(@(p)p.period,nsbr_pred(1).point);
omegas2_pred=arrayfun(@(p)p.period,nsbr_pred(2).point);
plot(getpars(nsbr(1),ind.mu1),omegas1,'Color',cm(1,:));
plot(getpars(nsbr(2),ind.mu1),omegas2,'Color',cm(2,:));
plot(getpars(nsbr_pred(1),ind.mu1),omegas1_pred,'.','Color',cm(1,:));
plot(getpars(nsbr_pred(2),ind.mu1),omegas2_pred,'.','Color',cm(2,:));
xlabel('$\mu_1$','Interpreter','LaTex');
ylabel('$\mu_2$','Interpreter','LaTex');
title('Compare computed and predicted periods')
%% save matlab workspace to file
save('vdpo.mat')