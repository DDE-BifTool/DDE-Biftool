%% test branch point detection
%
% $Id$
%%
%% Add paths and load sym package if octave is used
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities');
format compact
format short g
%% Set number of delays and parameter names
% This could also be loaded from mat file |'minimal_demo_parnames.mat'|.
parnames={'p','q','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
funcs=set_symfuncs(@sym_bp_demo,'sys_tau',@()ip.tau);
x0=[0;0.1]
par0([ip.p,ip.q,ip.tau])=[-1,x0(2)+x0(2)^3,1];
bds={'max_bound',[ip.tau,3;ip.p,1],'min_bound',[ip.p,-2]};
[stst,suc]=SetupStst(funcs,'x',x0,'parameter',par0,'contpar',ip.p,'dir',ip.p,...
    bds{:})
figure(1);clf;
stst=br_contn(funcs,stst,100);
nunst=GetStability(stst,'funcs',funcs);
%%
ibp=find(nunst==1,1,'first');
[bp,suc,blin]=br_bp_solve(funcs,stst,ibp,'print_residual_info',1)
%% branch off
stst2=SetupBranchSwitch(funcs,stst,ibp+(-1:0))
stst2=br_contn(funcs,stst2,100);
stst2=br_rvers(stst2);
stst2=br_contn(funcs,stst2,100);
%% periodic orbits
[stst2,suc]=SetupStst(funcs,'x',x0,'parameter',par0,'contpar',ip.tau,'dir',ip.tau,...
    bds{:})
figure(1);clf;
stst2=br_contn(funcs,stst2,100);
nunst2=GetStability(stst2,'funcs',funcs);
%%
ih=find(nunst2==2,1,'first');
[psol,suc]=SetupPsol(funcs,stst2,ih)
psol=br_contn(funcs,psol,100);
%%
[psolp,suc]=ChangeBranchParameters(funcs,psol,length(psol.point),'contpar',ip.p,'dir',ip.p);
psolp=br_contn(funcs,psolp,100);
nunstp=GetStability(psolp,'funcs',funcs,'exclude_trivial',true);
%%
ipbp=find(nunstp==1,1,'first');
[bpp,sucp,blinp]=br_bp_solve(funcs,psolp,ipbp,'print_residual_info',1)
[r,J]=p_correc_rhs(funcs,psolp.method.point,bpp,psolp.parameter.free);
assert(norm(r,'inf')<psolp.method.point.minimal_accuracy);
assert(rank(full(J))==size(J,1)-1);
assert(norm(J*dde_x_from_point(blinp.v0,ip.p))<psolp.method.point.minimal_accuracy);
assert(norm(blinp.w0'*J)<psolp.method.point.minimal_accuracy);
%%
psol2=SetupBranchSwitch(funcs,psolp,ipbp+(-1:0))
psol2=br_contn(funcs,psol2,30);
psol2=br_rvers(psol2);
psol2=br_contn(funcs,psol2,30);
