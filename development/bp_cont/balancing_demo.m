%%
clear
addpath([pwd(),'/../../ddebiftool']);
addpath([pwd(),'/../../ddebiftool_utilities'],...
    [pwd(),'/../../ddebiftool_extra_nmfm'],...
    [pwd(),'/../../ddebiftool_extra_psol'])    
parnames={'phi','psi','r','epsilon'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
funcs=set_symfuncs(@sym_balancing);
%% initial parameters
parini([ip.phi,ip.psi,ip.r,ip.epsilon])=...
      [  -0.2,    0.6, 0.5,   0];
xini=[0;0];
bounds={'max_bound',[ip.phi,1;ip.psi,1.1],...
    'min_bound',[ip.phi,-0.25;ip.psi,0.45]};
%% Follow equilibria, find Hopf bifurcation & its criticality
eqs=SetupStst(funcs,'x',xini,'parameter',parini,'step',0.1,...
    'contpar',ip.phi,'max_step',[ip.phi,0.01],bounds{:});
figure(1);clf;ax1=gca;xlabel(ax1,'phi');
eqs=br_contn(funcs,eqs,100,'plotaxis',ax1);
[nunst_eqs,dum,dum,eqs.point]=GetStability(eqs,'funcs',funcs);
%%
ibp=find(nunst_eqs==0,1,'first');
[bpfuncs,bpbr,suc]=SetupBPstst(funcs,eqs,ibp,[1,0],'contpar',[ip.phi,ip.psi],...
    'dir',ip.psi,'step',0.05,'print_residual_info',1);
figure(2);clf;ax2=gca;xlabel(ax2,'phi');ylabel(ax2,'psi');
bpbr=br_contn(bpfuncs,bpbr,100,'plotaxis',ax2);
bpbr=br_rvers(bpbr);
bpbr=br_contn(bpfuncs,bpbr,100,'plotaxis',ax2);
%%
[eqs_wbifs,testfuncs,bifind,biftype]=LocateSpecialPoints(funcs,eqs,...
    'debug',true,'print',1,'extra_cond',[1,0]);
%% Continue Hopf bifurcations
ihopf=find(strcmp(biftype,'hopf'));
hopf=SetupHopf(funcs,eqs_wbifs,bifind(ihopf),'contpar',[ip.phi,ip.psi],'dir',ip.psi,...
    bounds{:},'max_step',[0,2e-2]);
figure(2);
hopf=br_contn(funcs,hopf,50,'plotaxis',ax2);
hopf=br_rvers(hopf);
hopf=br_contn(funcs,hopf,30,'plotaxis',ax2);
%% Find degenerate Hopf bifurcations
[hopf_wbifs,hopffuncs,bif2ind,bif2type]=LocateSpecialPoints(funcs,hopf,'print',1);
%% plot 2d bif
figure(2);clf;ax2=gca;hold(ax2,'on');xlabel(ax2,'phi');ylabel(ax2,'psi');
Plot2dBranch(bpbr,'ax',ax2,'funcs',bpfuncs);
Plot2dBranch(hopf_wbifs,'ax',ax2);
%% branch off symmetric periodic orbit
% The branch will connect to the other degenerate Hopf bifurcation.
figure(1);
[perbranch,suc]=SetupPsol(funcs,eqs_wbifs,bifind(ihopf),...
    'degree',6,'intervals',50,'print_residual_info',1,'max_step',[0,1]);
%%
perbranch=br_contn(funcs,perbranch,100,'plotaxis',ax1);
%% Check Stability for POs
[nunst_per,dom_per,triv_per,perbranch.point]=GetStability(perbranch,'funcs',funcs,...
    'exclude_trivial',true);
%% find symmetry breaking
isb=find(diff(nunst_per)==1,1,'first')+1;
[sbfuncs,sbbranch,suc]=SetupPOfold(funcs,perbranch,isb,'nextrapar',1,...
    'extra_cond',{@(p,pref)sys_cond_POBP(p,[1,0])},...
    'contpar',[ip.phi,ip.psi],'dir',ip.phi,'step',0.1,'max_step',[],...
    'minimal_angle',0.6);
sbbranch=br_contn(sbfuncs,sbbranch,40,'plotaxis',ax2);
sbbranch=br_rvers(sbbranch);
sbbranch=br_contn(sbfuncs,sbbranch,700,'plotaxis',ax2);
%%
[nunst_sb,dom_sb,triv_sb,sbbranch.point]=GetStability(sbbranch,'funcs',sbfuncs,'exclude_trivial',true);
%% Plot complete bifurcation diagram so far
figure(2);clf;ax=gca;hold(ax,'on');xlabel(ax,'phi');ylabel(ax,'psi');
Plot2dBranch(hopf_wbifs,'ax',ax);
Plot2dBranch(bpbr,'ax',ax,'funcs',bpfuncs,'lgname','pitchfork');
Plot2dBranch(sbbranch,'ax',ax,'funcs',sbfuncs,'lgname','PO pitchfork');
set(ax,'fontsize',18,'fontweight','bold','fontname','courier','linewidth',2,'box','on');
%ax.YLim(2)=6.5;