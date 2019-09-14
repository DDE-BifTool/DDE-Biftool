%% Coupled FHN neural system with delay
% 
% Demo illustrating how to branch off a generalied Hopf bifurcation point
%
%% Differential equations coupled FHN neural system with delay
% 
% From: Zhen, Bin and Xu, Jian
% Bautin bifurcation analysis for synchronous solution of a coupled FHN
% neural system with delay.
% Commun. Nonlinear Sci. Numer. Simul. 15 (2010), no. 2, 442–458.
%
% $$u_1'=-\frac{u_1^3}{3}+(c+\alpha)u_1^2+d u_1-u_2+2\beta*\text{tanh}(u_1(t-\tau)),$$
%
% $$u_2'=\epsilon(u_1-b u_2).$$
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
% This could also be loaded from mat file |'FHN_parnames.mat'|.
parnames={'beta','alpha','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Set the funcs structure
% We load the precalculated multilinear forms. These have been generated
% with the file gen_sym_FHN.m.
funcs=set_symfuncs(@sym_FHN_mf,'sys_tau',@()ind.tau);
% Alternatively, uncomment the line below to use directional derivaties.
% funcs=set_symfuncs(@sym_FHN,'sys_tau',@()ind.tau);
%% Generalized-Hopf point derived in the article of Zhen, Bin and Xu, Jian
% construct steady-state point
beta0=1.9; alpha0=-0.9710; tau0=1.7722;
stst=dde_stst_create('x',[0;0]);
stst.parameter(ind.beta)  = beta0;
stst.parameter(ind.alpha) = alpha0;
stst.parameter(ind.tau)   = tau0;
% Calculate stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1(1:6)
%% Calculate critical normal form coefficients
hopf=p_tohopf(funcs,stst);
method=df_mthod(funcs,'hopf');
hopf.stability=p_stabil(funcs,hopf,method.stability);
genh=p_togenh(hopf);
genh=nmfm_genh(funcs,genh);
genh.nmfm.L1
% not a generalizd-Hopf point!
%% Initialize Hopf branch
unfolding_pars=[ind.beta, ind.alpha];
hbr=df_brnch(funcs,unfolding_pars,'hopf');
hbr.point=hopf;
hbr.point(2)=hopf;
hbr.point(2).parameter(ind.alpha)=...
hbr.point(2).parameter(ind.alpha)+0.001;
method=df_mthod(funcs,'hopf');
method.point.print_residual_info=1;
hbr.point(2)=p_correc(funcs,hbr.point(2),ind.beta,[],method.point);
%% Continue Hopf branch
figure(1); clf;
hbr=br_contn(funcs,hbr,30);
hbr=br_rvers(hbr);
hbr=br_contn(funcs,hbr,30);
title('Continued Hopf branch'); 
xlabel('$\beta$','Interpreter','LaTex')
ylabel('$\alpha$','Interpreter','LaTex')
box on
%% Criticality along Hopf bifurcation, codimension 2 points and their normal forms
% The function LocateSpecialPoints can be controlled and performed manually (see
% <minimal_demo_extra_nmfm.html>). The field |hopftests.genh| contains the
% criticality of the Hopf bifurcation.
[hbr_wbifs,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hbr);
al=arrayfun(@(x)x.parameter(ind.alpha),hbr_wbifs.point);
figure(2); clf;
plot(al,hopftests.genh(1,:),'.-',al,zeros(size(al)));
title('Criticality along Hopf bifurcation curve')
xlabel('$\beta$','Interpreter','LaTex');
ylabel('First Lyapunov coefficient (L1)')
%% Find branches of Codimension 1 bifurcations emerging from special points
% This helper function performs initialization only. By default, it finds
% the first 3 points along each branch.  The see
% <../../../ddebiftool_utilities/html/BranchFromCodimension2.html> for
% description of optional inputs and fields of output.
%  C1info=BranchFromCodimension2(funcs,hbr_wbifs);
%% Continue newly found branch
% npfold=50; lpcbr=br_contn(C1info(1).funcs,C1info(1).branch,npfold);
%% Calculate parameter-dependent normal form coefficients
% Select the located generalized Hopf point on the hopf_br_wbifs.
% Then convert the Hopf point to a generalized Hopf point.
% As before, we use the function nmfm_genh to calculate the normal form
% coefficients. By adding the option free_pars and providing the 
% unfolding parameters the parameter-dependent normal form coefficients
% are calculated.
hopf=hbr_wbifs.point(hc2_indices); % select generelized Hopf point
genh=p_togenh(hopf);
genh=nmfm_genh(funcs,genh,'free_pars',unfolding_pars);
genh.nmfm.L1
%% Continue Limit Point of Cycles curve emanating from generalized-Hopf point
figure(1); [lpcfuncs,lpcbr,~]=C1branch_from_C2point(funcs,genh,...
    unfolding_pars,'codim2',genh.kind,'codim1','POfold');
nop=50; [lpcbr,suc]=br_contn(lpcfuncs,lpcbr,nop); assert(suc>0)
%% Predictor Limit Point of Cycles curve emanating from generalized-Hopf point
[~,lpcbr_pred]=C1branch_from_C2point(funcs,genh,unfolding_pars,...
    'codim2',genh.kind,'codim1','POfold',...
    'step',linspace(0,1,45),'predictor',1);
%% Bifurcation diagram
figure(3); clf; hold on;
% Inline function to subtract parameters
getpars=@(points,ind) arrayfun(@(p)p.parameter(ind),points);
cm=colormap('lines');
L1s=hopftests.genh(1,:); % First Lyapunov coefficients along the hopf branch
% Plot sub- and supercritical Hopf branches
plot(getpars(hbr_wbifs.point(L1s>0),ind.beta),...
     getpars(hbr_wbifs.point(L1s>0),ind.alpha),'Color',cm(1,:),...
     'DisplayName','subcritical Hopf branch');
 plot(getpars(hbr_wbifs.point(L1s<0),ind.beta),...
     getpars(hbr_wbifs.point(L1s<0),ind.alpha),'Color',cm(2,:),...
     'DisplayName','supercritical Hopf branch');
plot(getpars(lpcbr.point,ind.beta),...
     getpars(lpcbr.point,ind.alpha),'Color',cm(3,:),...
     'DisplayName','LPC branch');
plot(getpars(lpcbr_pred.point,ind.beta),...
     getpars(lpcbr_pred.point,ind.alpha),'.','Color',cm(3,:),...
     'DisplayName','LPC predictor');
plot(getpars(genh,ind.beta),getpars(genh,ind.alpha),'k.',...
    'MarkerSize',8,'DisplayName','generalized Hopf point');
title('Bifurcation diagram near generalized Hopf point')
xlabel('$\beta$','Interpreter','LaTex');
ylabel('$\alpha$','Interpreter','LaTex');
text(1.8779,-1.1001,'I');
text(1.9130,-0.9008,'II');
text(1.8890,-0.5854,'III');
axis([1.86 1.945 -1.5000 -0.4000])
legend('Location','NorthEast');
box on;
%% Plot comparing computed and predicted periodic orbits
figure(4); clf; hold on;
for i=1:15
   plot(lpcbr.point(i).profile(1,:),...
       lpcbr.point(i).profile(2,:),'Color',cm(1,:));
end
% Plot predicted periodic orbits
for i=1:7
   plot(lpcbr_pred.point(i).profile(1,:),...
       lpcbr_pred.point(i).profile(2,:),'Color',cm(2,:));
end
% add labels
xlabel('$u_1$','Interpreter','LaTex');
ylabel('$u_2$','Interpreter','LaTex');
% add title
title('Compare computed and predicted periodic orbits')
box on
%% Compare computed and predicted periods
figure(5); clf; hold on;
% plot computed periods on nsbr(1) and nsbr(2)
omegas1=arrayfun(@(p)p.period,lpcbr.point);
omegas1_pred=arrayfun(@(p)p.period,lpcbr_pred.point);
plot(getpars(lpcbr.point,ind.alpha),omegas1,'Color',cm(1,:));
plot(getpars(lpcbr_pred.point,ind.alpha),omegas1_pred,'.','Color',cm(1,:));
title('Compare computed and predicted periods')
xlabel('$\alpha$','Interpreter','LaTex')
ylabel('$\omega$','Interpreter','LaTex')
axis([-1.1024   -0.3442   24.4409  201.7154])
legend('Computed period','Predicted period','Location','SouthEast'); box on
%% Save and continue to simulation near generalized Hopf point <simulation.html>
save('FHN_results.mat');