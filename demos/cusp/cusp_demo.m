%% Demo - System with cusp and Takens-Bogdanov bifurcations
%
% Demo contributed by Maikel Bosschaert, extended and modified by Jan
% Sieber
% 
% from Giannakopoulos, F. and Zapp, A. (2001). Bifurcations in a planar
% system of differential delay equations modeling neural activity. Physica
% D: Nonlinear Phenomena, 159(3):215-232.
%
%% Differential equations
%
% 
% $$\mu\dot{u}_{1}(t)=-u_{1}(t)+q_{11}\alpha(u_{1}(t-T))-q_{12}u_{2}(t-T)+e_{1}$$
%
% $$\mu\dot{u}_{2}(t)=-u_{2}(t)+q_{21}\alpha(u_{1}(t-T))-q_{22}u_{2}(t-T)+e_{2}$$
%
%% Demo walk-through
% The demo will show how one can perform tasks related to equilibria and
% their bifurcations at various levels of automation. Functions in the
% folder |ddebiftool_utilities| and the function |br_bifdet| provide
% automatic detection of bifurcations and computation of normal forms.
%
% However, these automated functions are not very robust such that the
% step-by-step approach is illustrated, too.
%
% <html>
% $Id: cusp_demo.m 176 2017-03-13 00:25:33Z jansieber $
% </html>
%
%% Definition of paths, problem right-hand side, delays and derivatives
% Note that derivatives for the problem have been defined manually in
% |cusp_mfderi|. Removing the name-value pair for |'sys_mfderi'| (or
% |'sys_deri'|) will make DDE-Biftool switch to finite difference
% approximations instead.

%#ok<*NOPTS>
clear      % clear variables
close all; % close figures
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');

disp('Cusp demo');
format compact
%% Illustrations of choices for problem definition
% There are several ways to specify a problem (the DDE) for DDE-Biftool
% when the delays are parameters. The simplest ways is to just specify the
% right-hand side in a function |sys_rhs(xx,par)|, where |xx(i,j)| is
% $x_i(t-\tau_j)$, and the function |sys_tau()| which returns the indices
% of the delays in the parameter vector. All derivatives will be
% approximated numerically in this case.
fminimal=set_funcs('sys_rhs',@cusp_rhs,'sys_tau',@()6,'x_vectorized',true);
%% Manual provision of derivatives
% Alternatively, one may set |'sys_mfderi'| to |@cusp_mfderi| to define
% high-order derivatives with respect to the state in equilibria or
% |'sys_deri'| to |@cusp_deri|. This is the way of specifying derivatives
% that has been supported by DDE-Biftool traditionally.
ftraditional=set_funcs('sys_rhs',@cusp_rhs,'sys_tau',@()6,...
    'sys_deri',@cusp_deri,'sys_mfderi',@cusp_mfderi);
%% Directional derivatives, manually provided
% A new way of providing derivatives manually is to provide the derivative
% in (xx,par) a given direction (dx,dp), where dx and q have the same
% format as xx and par. Derivatives up to order 2 are used for standard
% bifurcation analysis (continuation of steady states, periodic orbits and
% codimension-one bifurcations) and derivatives up to order 5 (only with
% dp=0) are useful for normal form computations. If only lower-order
% derivatives are provided, they will be called also for finite difference
% approximations of higher-order derivatives. A recommended format is a
% cell array |dirderi|, where |dirderi{k}| is the kth-order derivative. For
% this example, only the first two orders were hand-coded in the function
% file |cusp_dirderi|:
dirderi{1}=@(xx,par,dx,dp)cusp_dirderi(1,xx,par,dx,dp);
dirderi{2}=@(xx,par,dx,dp)cusp_dirderi(2,xx,par,dx,dp);
fdirectional=set_funcs('sys_rhs',@cusp_rhs,'sys_tau',@()6,...
    'sys_dirderi',dirderi,'x_vectorized',true);
%% Symbolic Toolbox
% The script <gen_sym_cusp.html> shows how one can use the symbolic toolbox
% to generate all derviatives in a function file |sym_cusp_rhs|. Due to
% limitations and quirks of the symbolic toolbox (as of version 7.0) some
% wrappers around the generated functions are needed. The wrapper
% |set_symfuncs| around |set_funcs| provides them. The delays still have to
% be specified.
fsymbolic=set_symfuncs(@sym_cusp,'sys_tau',@()6);
%% Select which problem definition to test
funcs=fsymbolic;
%% Continuation of equilibria in the parameter e1 (E)
% The order of parameters in the vector |par| in the right-hand side is
% determined by the assignment to |par|.
Q=0;
E=0;

q11=2.6; q12=Q; q21=1; e1=E; e2=0; T=1;
par = [q11,q12,q21,e1,e2,T];

stst_br = SetupStst(funcs,'x',[0;0],'parameter',par,...
    'contpar',4,'max_step',[4 0.02],'max_bound',[4 1],'min_bound',[4 0],...
    'newheuristics_tests',0);
nsteps=300;
figure(1);clf;ax1=gca;
stst_br = br_contn(funcs,stst_br,nsteps,'plotaxis',ax1)
disp('stst_br.point = '); disp(stst_br.point);
%% Manual finding of stability
% After continuation the field |'stability'| of the entries in
% |stst_br.point| is still empty. We can obtain
% stability for each point manually using the core routine |p_stabil|. Here
% it is done for the first element of stst_br.point. The parameters for
% stability computation are set in the structure |stst_br.method.stability|
% (see manual).
disp('stst_br.point(1).stability = '); disp(stst_br.point(1).stability);
stability=p_stabil(funcs,stst_br.point(1),stst_br.method.stability)
%% Stability along branch
% Determining the stability along an entire branch can be done with either
% the standard function |br_stabl| or with the convenience fucntion
% GetStability. GetStability also immediately extracts the number of
% unstable eigenvalues and the dominant eigenvalue.
stst_br=br_stabl(funcs,stst_br,0,0);
[stst_nunst,stst_dom,~,stst_br.point]=GetStability(stst_br,'funcs',funcs);
disp('number of unstable eigenvalues = '); disp(stst_nunst');
%% Automated finding of stability, bifurcations and their normal forms
% The highest level of automatino is provided by the convenience function
% |LocateSpecialPoints|. It determines stability (if not yet done), detects
% bifurcations and computes their normal forms. The outputs are the branch
% |stst_br| where now all fields have complete stability information and
% where bifurcation points have non-empty flag fields together with nvec
% and nmfm information. The field |nmfm| contains the normal form
% information. |stst_tests| is a stucture containing a field for each test
% function applied along the branch. |codim1_ind| are the indices in the
% |stst_br.point| array at which bifurcations were encountered. |stst_bif|
% is a cell array of strings of the same length as indices, containing the
% type of the bifurcations.
[stst_br,stst_tests,codim1_ind,stst_bif]=LocateSpecialPoints(funcs,stst_br) 
%% Displaying branch data
% Parameters and solution measures are available through the core functions
% br_measr and df_measr, or they can be extracted manually:
E=arrayfun(@(p)p.parameter(stst_br.parameter.free(1)),stst_br.point);
xval=arrayfun(@(p)p.x(1),stst_br.point);
cla(ax1);
plot(ax1,E(stst_nunst==0),xval(stst_nunst==0),'b.',...
    E(stst_nunst>0),xval(stst_nunst>0),'r.');
xlabel(ax1,'E');
ylabel(ax1,'x1');
%% Convenience plot
% Alternatively, the function |Plot2dBranch| plots a branch with various default setting
% that can be overwritten by optional name-value pairs.
Plot2dBranch(stst_br,'ax',ax1);
xlabel(ax1,'E');
ylabel(ax1,'x1');
%% Find fold and continue it in (E,Q)
% The fold detected (for example, by LocateSpecialPoints) can be tracked in
% two parameters using the convenience function |SetupFold| and then the
% core function |br_contn|.
Qmin=-0.1; Qmax=2; Emin=-0.6; Emax=0.6;
bifbounds={'min_bound',[2 Qmin; 4 Emin],'max_bound',[2 Qmax; 4 Emax],...
    'max_step',[2 0.02; 4 0.02]};
[fold_branch, suc] = SetupFold(funcs, stst_br, codim1_ind(1), ...
    'contpar', [2 4], 'dir', 4, 'step', 0.02,bifbounds{:});


figure(2);clf;ax2=gca;
fold_branch=br_contn(funcs,fold_branch,300,'plotaxis',ax2);
fold_branch = br_rvers(fold_branch);
fold_branch=br_contn(funcs,fold_branch,300,'plotaxis',ax2);
%% Manual detection and computation of codimension two bifurcations
% The option |'exclude_trivial'| permits one to ignore the trivial
% eigenvalue along a fold branch such that stability changes transversal to
% the known center direction do not get obscured. The output |triverr|
% shows the difference between the eigenvalues closest to zero and the
% theoretically known value (0).
%
% Changes by 1 suggest Bogdanov-Takens bifurcations, changes by 2 suggest
% Zero-Hopf (Gavrilov-Guckenheimer) bifurcations.
[nunst_fold,dom_fold,triverr,fold_branch.point]=GetStability(fold_branch,...
    'funcs',funcs,'exclude_trivial',true);
disp('number of transversally unstable eigenvalues of fold = ');
disp(nunst_fold');
%% Manual computation of Bogdanov-Takens (BT) bifurcations
% There are two changes indicating BT. The convenience function
% TakensBogdanovNormalform refines the branch and computes the normal form.
% The refined branch is its 3rd output (not shown here).
ind_BT=find(abs(diff(nunst_fold))==1)
bt(1)=TakensBogdanovNormalform(funcs,fold_branch,ind_BT(1)+(0:1));
bt(2)=TakensBogdanovNormalform(funcs,fold_branch,ind_BT(2)+(0:1));
celldisp(...
    arrayfun(@(p)p.parameter(fold_branch.parameter.free),bt,'uniformoutput',false),...
    'parameters at BT points ')
celldisp({bt.nmfm},'normal form coefficients at BT points ')
%% Computation of fold normal form coefficients along branch
% Sign changes of the fold normal form coefficient a2 can be used to detect
% a cusp (which is also visible in the online plot as the sharp  point at
% maximum value of parameter 2). Note that a2 has a singularity at BT
% points.
a2=FoldNormalformCoefficients(funcs,fold_branch);
E_fold=arrayfun(@(p)p.parameter(fold_branch.parameter.free(2)),fold_branch.point);
Q_fold=arrayfun(@(p)p.parameter(fold_branch.parameter.free(1)),fold_branch.point);
figure(3);
plot(E_fold,a2,'.-');
set(gca,'ylim',[-3,3]);
grid on
%% Manual computation of cusp
% The function CuspNormalform refines the branch and computes the normal
% form coefficient at the cusp.
ind_cu=find(diff(sign(a2)));
disp('sign change of a2 at index = '); disp(ind_cu);
cu=CuspNormalform(funcs,fold_branch,ind_cu(2)+(0:1));
disp('parameters at cusp = '); disp(cu.parameter(fold_branch.parameter.free));
disp('normal form coefficients at cusp = '); disp(cu.nmfm);
%% Automatically detect and compute codimension 2 bifurcations: cusp and Bogdanov-Takens
% The above tasks are automated using the convenience function
% |LocateSpecialPoints|. The first output is the refined branch where
% special poitns are indicated by non-empty flags and all points have
% either their fold normal form coefficient or their codimension 2 normal
% form. The second output is a struct with all relevant test function
% values along the branch. The 3rd output contains the indices in
% |fold_branch.point| at which special points were detected. The subsequent
% call to GetStability extracts stability information and computes it anew
% where necessary. The count of unstable eigenvalues takes into account the
% number of known eigenvalues at BT points.
[fold_branch,fold_tests,codim2_ind,codim2types]=LocateSpecialPoints(funcs,fold_branch);
fold_nunst=GetStability(fold_branch,'exclude_trivial',true);
%% Automated plot of fold branch in two parameters
% The function Plot2dBranch plots the branch, split into regions of
% different stability and highlighting spcial points with reasonable
% default values. Its first output are the legend references and entries,
% which can be passed on as input when adding branches to the plot.
lg2=Plot2dBranch(fold_branch,'ax',ax2);
grid on
xlabel(ax2,'Q');
ylabel(ax2,'E');
%% Compute additional line of equilibria at E=0
% This branch passes through an apparent fold, which is, however, a
% pitchfork at E=0 (the cusp in the (E-Q) plane.
stst2 = ChangeBranchParameters(funcs,stst_br,length(stst_br.point),'contpar',2,...
    'max_step',[0 0.05],'max_bound',[2 2],'min_bound',[2 0],'step',0.01);
nsteps=300;
figure(3);clf;ax3=gca;
stst2 = br_contn(funcs,stst2,nsteps,'plotaxis',ax3);
%% Detect bifurcations along branch and determine stability
% Note that the Hopf bifurcation is close to a BT point such that the
% Newton iteration for its defining system may take longer to converge. We
% monitor this by setting |print_residual_info|.
[stst2,~,hopf_ind,stst2_types]=LocateSpecialPoints(funcs,stst2,'print_residual_info',1);
%% Continue Hopf bifurcations
% The system is symmetric such that both Hopf bifurcation curves are mirror
% images of each other in the parameter plane. We continue only one of them.
figure(4);clf;ax4=gca;
for i=length(hopf_ind):-1:1
    hopf{i}=SetupHopf(funcs,stst2,hopf_ind(i),'contpar',[2,4], 'dir', 4, 'step', 0.02,...
        'print_residual_info',0,bifbounds{:});
    hopf{i}=br_contn(funcs,hopf{i},30,'plotaxis',ax4);
    hopf{i}=br_rvers(hopf{i});
    hopf{i}=br_contn(funcs,hopf{i},30,'plotaxis',ax4);
    [hopf{i},hopf_tests{i},hopfc2_ind{i},hopfc2types{i}]=LocateSpecialPoints(funcs,hopf{i});
end
%% Plot resulting 2-parameter bifurcation diagram for equilibria
% Note that the Hopf bifurcations are subcritical close to the BT points.
hold(ax2,'on');
lg2=Plot2dBranch(hopf{1},'ax',ax2,'oldlegend',lg2);
lg2=Plot2dBranch(hopf{2},'ax',ax2,'oldlegend',lg2);
set(ax2,'xlim',[1,2]);
legend(lg2{1},lg2{2},'location','east');
%% Save results
save('cusp_demo_results.mat');