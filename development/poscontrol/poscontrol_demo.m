%%
clear
ntau=3;
parnames={'tau0','s0','K','c'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%rhs=@(x,s,tau0,s0,c,K)[-K*c/2*(s(2,:)-s0);...
%    (-K*c/2*(s(4,:)+s(2,:)-2*s0)-c*s(1,:)+x(3,:)+x(1,:))/c+K*c/2*(s(4,:)-s0)];
%%rhs_xp=@(xx,par)rhs(reshape(xx(1,:,:),ntau+1,[]),reshape(xx(2,:,:),ntau+1,[]),par(1),par(2),par(3),par(4));
%funcs=set_funcs('sys_rhs',@poscontrol_rhs,'sys_tau',@poscontrol_tau,'sys_ntau',@()3,...
%    'x_vectorized',true);
%%
funcs=set_symfuncs(@sym_poscontrol);
%% initial parameters
tau0=0.1;
s0=1;
c=1;
K=2;
x0=[1;1];
par0=cellfun(@(x)evalin('caller',x),parnames); %[beta,n,tau,gamma];
contpar=ind.tau0;
%%
eqs=SetupStst(funcs,'x',x0,'parameter',par0,'step',0.1,...
    'contpar',contpar,'max_step',[contpar,0.05],'max_bound',[contpar,4]);
figure(1);clf;ax1=gca;
eqs=br_contn(funcs,eqs,100,'plotaxis',ax1);

[eqs,testfuncs,bifind,biftype]=LocateSpecialPoints(funcs,eqs);
%%
%% branch off periodic orbits
psol0=SetupPsol(funcs,eqs,bifind,'contpar',ind.tau0,...
    'print_residual_info',1);
psol0=br_contn(funcs,psol0,20,'plotaxis',ax1);
[nunst0,dom0,triv0,psol0.point]=GetStability(psol0,'funcs',funcs,'exclude_trivial',true);

%%
hopf=SetupHopf(funcs,eqs,bifind(1),'contpar',[ind.tau0,ind.s0],'dir',ind.s0,...
    'min_bound',[ind.tau0,0;ind.s0,0],'step',1e-3,'max_step',[0,1e-1]);
figure(2);clf;ax2=gca;
hopf=br_contn(funcs,hopf,100,'plotaxis',ax2);
hopf=br_rvers(hopf);
hopf=br_contn(funcs,hopf,100,'plotaxis',ax2);
%%
[hopf_wbifs,hopffuncs,bif2ind,bif2type]=LocateSpecialPoints(funcs,hopf);
%% plot 2d bif
figure(2);clf;
[lg2,ax2]=Plot2dBranch(hopf_wbifs);
%% branch off periodic orbits
psol1=SetupPsol(funcs,hopf_wbifs,bif2ind-20,'contpar',ind.tau0);
psol1=br_contn(funcs,psol1,20,'plotaxis',ax1);
[nunst1,dom1,triv1,psol1.point]=GetStability(psol1,'funcs',funcs,'exclude_trivial',true);
psol2=SetupPsol(funcs,hopf_wbifs,bif2ind+20,'contpar',ind.tau0);
psol2=br_contn(funcs,psol2,20,'plotaxis',ax1);
[nunst2,dom2,triv2,psol2.point]=GetStability(psol2,'funcs',funcs,'exclude_trivial',true);
%% continue fold of p.o.
ifold=find(abs(diff(nunst1))==1,1,'last');
[pfuncs,pfold]=SetupPOfold(funcs,psol1,ifold,'print_residual_info',1,...
    'contpar',[ind.tau0,ind.s0],'step',-1e-2,'dir',ind.s0);
pfold=br_contn(pfuncs,pfold,20,'plotaxis',ax2);
pfold=br_rvers(pfold);
pfold=br_contn(pfuncs,pfold,20,'plotaxis',ax2);
pfold=br_contn(pfuncs,pfold,20,'plotaxis',ax2);
[nunst_pf,dom_pf,triv_pf,pfold.point]=GetStability(pfold,'funcs',pfuncs,'exclude_trivial',true);
%%
figure(2);clf;ax2=gca;hold(ax2,'on');
[lg2,ax2]=Plot2dBranch(hopf_wbifs,'ax',ax2);
[lg2,ax2]=Plot2dBranch(pfold,'ax',ax2,'oldlegend',lg2,'funcs',pfuncs);
set(ax2,'fontweight','bold');
%grid(ax2,'on');
axis(ax2,'tight');
xlabel('tau');
ylabel('sref');
%%
figure(3);ax3=gca;
tauhopf=arrayfun(@(x)x.parameter(ind.tau0),hopf_wbifs.point);
L1hopf=hopffuncs.genh;
plot(tauhopf,L1hopf,'.-',[1,max(tauhopf)],[0,0],'k-',tauhopf(bif2ind),0,'ko',...
    'markerfacecolor','k');
set(ax3,'fontweight','bold');
xlim(ax3,[min(tauhopf),max(tauhopf)]);
xlabel('tau');
ylabel('L1');
%%
%save('poscontrol_results.mat');