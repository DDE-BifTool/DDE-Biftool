%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: tutorial_V_part2.m 134 2016-09-12 11:10:44Z mmbosschaert $
%% Tutorial V part II : Continuation of homoclinic orbits
%% Define system
clear variables
close all
addpath(...
    '../../../ddebiftool',...
    '../../../ddebiftool_utilities',...
    '../../../ddebiftool_extra_psol',...
    '../../../ddebiftool_extra_nmfm',...
    '../IV','../III');
% load results from tutorial III
load('../III/tutorial_III.mat');
HollingTanner_rhs = @(xx,par) [
        (xx(1,1,:)+par(4)).*(1-xx(1,1,:)-par(4))-...
            xx(1,1,:).*xx(2,1,:)./(par(3)*xx(2,1,:)+xx(1,1,:))-par(5);
        par(6)*xx(2,1,:).*(par(1)-xx(2,2,:)./xx(1,2,:))];
tau_ind=2;
funcs=set_funcs(...
    'sys_rhs', HollingTanner_rhs,...
    'sys_tau', @() tau_ind,...
    'sys_deri',@HollingTanner_deri,...
    'sys_mfderi', @HollingTanner_mfderi,...
    'x_vectorized',true);
%% Branch off to periodic orbit at some Hopf point, continue to large period
psol_branch=SetupPsol(funcs,hopf_branch_wbifs,2,...
    'contpar',inddelta,'degree',3,'intervals',50,parameter_bd{1:4},'max_step',[0,inf]);
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
%% plot period orbits
cm=colormap('lines');
pmeshes=cell2mat(arrayfun(@(x)x.mesh(:),psol_branch.point,'uniformoutput',false));
pprofs1=cell2mat(arrayfun(@(x)x.profile(1,:)',psol_branch.point,'uniformoutput',false));
pprofs2=cell2mat(arrayfun(@(x)x.profile(2,:)',psol_branch.point,'uniformoutput',false));
figure(5);clf
hold on
for i=1:length(pprofs1(1,:))
    plot(pprofs1(:,i),pprofs2(:,i),'color',cm(1,:));
end
xlabel('$x$','Interpreter','LaTex');
ylabel('$y$','Interpreter','LaTex');
plot(hopf_branch_wbifs.point(2).x(1),hopf_branch_wbifs.point(2).x(2),'.r')
drawnow;
%% Convert point close to end of the psol_branch to homoclinic orbit
% and correct homoclinic orbit.Then refine and repeat correction.
hcli=p_tohcli(funcs,psol_branch.point(end-5));
figure(6);clf;
p_pplot(hcli);
xlabel('time/period')
ylabel('$x,y$','Interpreter','LaTex');
drawnow;
mhcli=df_mthod(funcs,'hcli');
[hcli,suc]=p_correc(funcs,hcli,inddelta,[],mhcli.point); % correct
disp(suc);
hcli=p_remesh(hcli,3,50); % remesh it and
[hcli,suc]=p_correc(funcs,hcli,indbeta,[],mhcli.point); % correct it again
disp(suc);
%% Continue branch of homoclinic orbits in two parameters
hcli_br=df_brnch(funcs,[inddelta, indbeta],'hcli');
hcli_br.point=hcli;
hcli2=hcli;
hcli2.parameter(indbeta)=hcli2.parameter(indbeta)-1e-4;
[hcli2,suc]=p_correc(funcs,hcli2,inddelta,[],mhcli.point);
hcli_br.point(2)=hcli2;
hcli_br.parameter.max_bound=fold_branch_wbifs.parameter.max_bound;
hcli_br.parameter.min_bound=fold_branch_wbifs.parameter.min_bound;
hcli_br.parameter.max_step=[indbeta,5e-3;inddelta,5e-3];
hcli_br.method.point.print_residual_info=0;
figure(7);
hcli_br=br_contn(funcs,hcli_br,32);
hcli_br=br_rvers(hcli_br);
hcli_br=br_contn(funcs,hcli_br,40);
%% Add homoclinic orbit to two-parameter bifurcation diagram
beta_fold=getpar(fold_branch_wbifs,indbeta);
delta_fold=getpar(fold_branch_wbifs,inddelta);
beta_hopf=getpar(hopf_branch_wbifs,indbeta);
delta_hopf=getpar(hopf_branch_wbifs,inddelta);
delta_bt1=bgetpar(hopf_branch_wbifs,inddelta,'BT');
beta_bt1=bgetpar(hopf_branch_wbifs,indbeta,'BT');
figure(8); clf;
plot(beta_fold,delta_fold,'k');
hold on;
plot(beta_hopf,delta_hopf,'color',cm(1,:));
plot(beta_bt1,delta_bt1,'r.','MarkerSize',8);
plot(getpar(hcli_br,indbeta),getpar(hcli_br,inddelta),'color',cm(2,:));
xlabel('$\beta$','Interpreter','LaTex');
ylabel('$\delta$','Interpreter','LaTex');
legend('fold','Hopf','BT','homoclinic','Location','NorthWest')
axis([0.47   0.51   0.4   0.7])