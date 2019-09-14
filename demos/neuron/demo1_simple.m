%% DDE-BIFTOOL demo 1 - Neuron, simplification using utilities
%
% <html>
% $Id: demo1_simple.m 374 2019-09-14 14:02:58Z jansieber $
% </html>
%
% This demo repeats the illustrative example using the convenience
% functions in |ddebiftool_utilities|. The system of
% delay differential equations is [Shay,99] again
% 
% $$\left\{\begin{array}{l}\dot{x_1}(t)=
% -\kappa x_1(t)+\beta \tanh(x_1(t-\tau_s))+a_{12}\tanh(x_2(t-\tau_2)) \\ 
% \dot{x_2}(t)=-\kappa x_2(t)+\beta \tanh(x_2(t-\tau_s))+a_{21}\tanh(x_1(t-\tau_1)).
% \end{array}\right.$$
%
% This system models two coupled neurons with time delayed connections.
% It has two components ($x_1$ and $x_2$), three delays 
% ($\tau_1$, $\tau_2$ and $\tau_s$), and four other parameters 
% ($\kappa$, $\beta$, $a_{12}$ and $a_{21}$).
%
% The demo will show (with fewer explanations than <demo1_funcs.html> and
% its follow-on files)
% 
% * continuation of equilibria in $a_{21}$,
% * computation of their stability (eigenvalues of linearization),
% * continuation of Hopf bifurcations in $a_{21}$ and $\tau_s$,
% * branching off from a Hopf bifurcation to continue a family of periodic
% orbits in $a_{21}$,
% * computation of stability of periodic orbits (Floquet multipliers of
% linearization),
% * continuation of folds of periodic orbits in $a_{21}$ and $\tau_s$
%%

%% Additional Folder for Path
% The utility functions are stored in a separate folder, which has to be
% loaded in addition to |ddebiftool|. We also include the folders for
% periodic orbits and for normal form functions.
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...
    [base,'ddebiftool_extra_psol/'],...
    [base,'ddebiftool_utilities/'],...
    [base,'ddebiftool_extra_nmfm/']);
clear;
format compact
format short g
%#ok<*ASGLU,*NOPTS,*NASGU>
%% Definition of right-hand side
% When using we vectorize the right-hand side in xx, but do not specify the
% derivatives analytically, relying on finite difference approximations.
% See <demo1_funcs.html> for details. The parameter vector has the order
% $[\kappa, \beta, a_{12}, a_{21}, \tau_1, \tau_2, \tau_s]$.
%
% When using fsymbolic we use the right-hand side and derivatives gnereated
% with the symbolic toolbox.
parnames={'kappa','beta','a12','a21','tau1','tau2','taus'};
cind=[parnames;num2cell(1:length(parnames))];
indp=struct(cind{:});
getp=@(br,name)arrayfun(@(x)x.parameter(indp.(name)), br.point);
neuron_sys_rhs=@(x,p)[...
    -p(indp.kappa)*x(1,1,:)+p(indp.beta)*tanh(x(1,4,:))+p(indp.a12)*tanh(x(2,3,:));....
    -p(indp.kappa)*x(2,1,:)+p(indp.beta)*tanh(x(2,4,:))+p(indp.a21)*tanh(x(1,2,:))];
neuron_tau=@()[indp.tau1,indp.tau2,indp.taus];
fnum=set_funcs(...
    'sys_rhs',neuron_sys_rhs,'sys_tau',neuron_tau,'x_vectorized',true);
fsymbolic=set_symfuncs(@sym_neuron,'sys_tau',neuron_tau);
%% Choose problem definition to be tested
funcs=fsymbolic;
% general continuation parameters, kept in a cell list
parbd={'min_bound',[indp.a21,0],'max_bound',[indp.a21,3; indp.taus,10],...
    'max_step',[indp.a21,0.2; indp.taus,0.5]};
%% Steady state branches
% The convenience function |SetupStst| can be used to define the initial
% piece of a steady-state branch. Its first arguments is |funcs|, the
% others are name-value pairs. Important parameters:
% 
% * |'parameter'|: row vector of initial parameters
% * |'x'|: column vector of initial equilibrium
% * |'contpar'|: index of continuation parameter (or vector of indices)
% * |'step'|: initial step along branch (default 0.01)
%
% All other name-value pairs are appended as fields to the structures in
% branch1 if their names match. Of course, the |branch| structure can also
% be manipulated manually afterwards. The subsequent continuation extends
% the branch in both directions up to the boundaries |min_bound| and
% |max_bound|.
par0([indp.kappa,indp.beta,indp.a12,indp.a21,indp.tau1,indp.tau2,indp.taus])=...
     [    1/2,      -1,      1,  2.34,       0.2      0.2,     1.5];
[branch1,suc]=SetupStst(funcs,...
    'parameter',par0,'x',[0;0],...
    'contpar',indp.a21,'step',0.1,parbd{:});
branch1.method.continuation.plot=0; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,100);
branch1=br_rvers(branch1);
[branch1,s,f,r]=br_contn(funcs,branch1,100)

%% Stability along branch and bifurcations
% The convenience function |GetStability| recomputes the eigenvalues if not
% yet present and returns as its first output |nunst| the number of
% unstable eigenvalues for bifurcation detection. Its first argument is the
% |branch| structure for which stability information is required. We find
% the first point at which the number of unstable eigenvalues changes by 2.
% The function |LocateSpecialPoints| also refines the branch as necessaryto
% pinpoint bifurcation points. Note that the branch point at a21=2.25 is
% incorrectly identified asa fold such that the correction fails.
branch1.method.stability.minimal_real_part=-2;
[branch1,bif1testfuncs,indbif1,bif1types]=LocateSpecialPoints(funcs,branch1);

%% Hopf continuation
% Similar to |SetupStst| the convenience function |SetupHopf| creates
% the initial Hopf branch. Its first arguments are |funcs|, the branch
% along which the Hopf bifurcation was detected (here |branch1|), and the
% index of the point near which the Hopf bifurcation was detected.
% Important parameters:
%
% * |'contpar'|: bifurcation parameters (vector of length >=2)
% * |'dir'|: index of parameter, which is varied at initial step. The
% default is [], which means that only one point on the branch is computed.
% This is useful if one wants to correct only a single Hopf point.
% * |'step'|: initial step along branch (default |1e-3|)
% * |'excudefreqs'|: list of frequencies that should be excluded (default
% []). The initial guess for the Hopf frequency is the complex conjugate
% pair closest to the imaginary axis, after one takes away a pair of
% eigenvalues for each frequency listed in |excludefreqs|.
%
% All other name-value pairs can be used to replace fields in the
% structures of the Hopf branch. Otherwise, the output |branch2| inherits
% all values from the input |branch|. The subsequent continuation  extends
% the branch in both directions up to the boundaries |min_bound| and
% |max_bound|.
indhopf=indbif1(strcmp(bif1types,'hopf'))
[branch2,suc]=SetupHopf(funcs,branch1,indhopf,'contpar',[indp.a21,indp.taus],...
    'dir',indp.taus,'step',0.1,'plot',1,parbd{:});
figure(1); clf;ax1=gca;
xlabel(ax1,'a21');ylabel(ax1,'\tau_s');
[branch2,s,f,r]=br_contn(funcs,branch2,40,'plotaxis',ax1);
branch2=br_rvers(branch2);
[branch2,s,f,r]=br_contn(funcs,branch2,20,'plotaxis',ax1);
%% Figure: Continuation (predictions and corrections) of Hopf bifurcation
% Predictions and corrections in the $(a_{21},\tau_s)$-plane after
% computation of a first branch of Hopf bifurcations.

%% Linear stability along Hopf curve
% The function |GetStability| has a few additional outputs and inputs. The
% optional input |exclude_trivial| forces exclusion of trivial eigenvalues
% for bifurcations (for example the pair on the imaginary axis for the Hopf
% bifurcation) in the count of unstable eigenvalues. The second output
% (here |dom|) contains the dominant eigenvalue (that is, the closest to
% the imaginary axis, excluding the trivial eigenvalues). The third output
% (here |triv_defect|) is non-empty whenever |exclude_trivial| is true. It
% contains the defect between the known trivial eigenvalue (say.
% 1i*point.omega for Hopf points) and the value obtained in the stability
% calculation. The fourth output is the modified |point| array, now
% containing stability information.
[hopf_branch_wbifs,hopftestfuncs,indbifhopf,bifhopftypes]=LocateSpecialPoints(funcs,branch2);
nunst_hopf=GetStability(hopf_branch_wbifs,'exclude_trivial',true);
% identify degenerate Takens-Bogdanov point (Hopf frequency crosses zero):
indTakens=indbifhopf(strcmp(bifhopftypes,'BT'));
pt_Takens=hopf_branch_wbifs.point(indTakens);
a21_Takens=pt_Takens.parameter(indp.a21);
%% Switch to second Hopf curve near double Hopf point
% Repeat the Hopf continuation for the point on |branch2| where another
% eigenvalue becomes unstable, again using |SetupHopf|, but now using
% $a_{21}$ as the parameter and initial step -0.05. Initially, we turn on
% residual printouts.
indhopf2=indbifhopf(find(strcmp('hoho',bifhopftypes),1,'first'));
[hopf2_branch,suc]=SetupHopf(funcs,hopf_branch_wbifs,indhopf2,'contpar',[indp.a21,indp.taus],...
    'dir',indp.a21,'step',-0.05,'print_residual_info',1);
hopf2_branch.method.point.print_residual_info=0;
[hopf2_branch,s,f,r]=br_contn(funcs,hopf2_branch,100,'plotaxis',ax1);
hopf2_branch=br_rvers(hopf2_branch);
[hopf2_branch,s,f,r]=br_contn(funcs,hopf2_branch,100,'plotaxis',ax1);
xlabel('a21');ylabel('\tau_s');
%% Figure: Continuation (predictions and corrections) of both Hopf bifurcations
% Predictions and corrections in the $(a_{21},\tau_s)$-plane after
% computation of second branch of Hopf bifurcations (superimposed on result
% of first Hopf bifurcation).
%% First Lyapunov coefficient L1 along both Hopf curves
% After a Hopf curve is computed, |HopfLyapunovCoefficients| returns the
% first Lyapunov coefficient L1, which determines if the Hopf bifurcation
% is supercritical (L1<0) or subcritical (L1>0). Since L1 depends on
% 3rd-order derivatives of |sys_rhs|, the default finite-difference
% approximation |df_mfderiv| is susceptible to round-off errors. Hence,
% |HopfLyapunovCoefficients| returns 2 results, the second one using a
% lower-order approximation. The difference between both outputs gives an
% error estimate. The results show that both curves are supercritical. When
% a symbolic derivative is provided the reported error will be zero (as it
% is not applicable).
[hopf2_branch_wbifs,hopf2testfuncs,indbifhopf2,bifhopf2types]=...
    LocateSpecialPoints(funcs,hopf2_branch);
fprintf('maximal L1 along branch2: %g\n',max(hopftestfuncs.genh(1,:)))
fprintf('  max error of L1 along branch2: %g\n',max(abs(diff(hopftestfuncs.genh))));
fprintf('maximal L1 along branch3: %g\n',max(hopf2testfuncs.genh(1,:)))
fprintf('  max error of L1 along branch3: %g\n',max(abs(diff(hopf2testfuncs.genh))));
%% Error in normal form of double Hopf point
% Similar to the L1 coefficients, the normal form coefficients of a
% Hopf-Hopf point depend on 3rd-order derivatives of |sys_rhs|. Thus, the
% function computing the normal form returns two values, the difference of
% which is an error estimate. The output of |HopfHopfNormalform| is an
% entire point structure of kind 'hoho'. The field nmfm contains normal
% form coefficients. Other fields of interest are |omega1| and |omega2|.
[hoho,hoho_low]=HopfHopfNormalform(funcs,hopf_branch_wbifs,indhopf2+[-1,1]);
fprintf(['Normal form coefficients of Hopf-Hopf point\n',...
    'at (a_21,tau_s)=(%g,%g) with omega1=%g, omega2=%g:\n'],...
    hoho.parameter(indp.a21),hoho.parameter(indp.taus),...
    hoho.nvec.omega(1),hoho.nvec.omega(2));
disp(hoho.nmfm);
fprintf('Error of normal form coefficients: %g\n',...
     norm(structfun(@(x)x,hoho.nmfm)-structfun(@(x)x,hoho_low.nmfm),'inf'));
%% Two-parameter bifurcation diagram of Hopf bifurcations
% with codimension-two points
figure(3);clf;
subplot(3,3,[1,2,4,5]);
ax3=gca;hold(ax3,'on');
args={'ax',ax3,'stability',0.75};
plot(ax3,a21_Takens*[1,1],[0,10],'r','linewidth',2,...
    'DisplayName','branch point','UserData','Plot2dBranch');
Plot2dBranch(hopf_branch_wbifs,args{:});
Plot2dBranch(hopf2_branch_wbifs,args{:});
set(findobj(get(ax3,'parent'),'type','legend'),'position',[0.7,0.1,0.2,0.2])
grid on
xlabel('a_{21}');
ylabel('\tau_s');
title('Two-parameter bifurcation diagram');
drawnow;
a21lim=get(ax3,'xlim');
tauslim=get(ax3,'ylim');
%% PLot Lyapunov coefficients
a21_h1 =getp(hopf_branch_wbifs,'a21');
taus_h1=getp(hopf_branch_wbifs,'taus');
a21_h2 =getp(hopf2_branch_wbifs,'a21');
taus_h2=getp(hopf2_branch_wbifs,'taus');
figure(3);
lw={'+-','linewidth',2};
subplot(3,3,[3,6]);
plot(hopftestfuncs.genh(1,:),taus_h1,lw{:});
set(gca,'ylim',tauslim,'xlim',[-0.6,0]);
grid on
ylabel('\tau_s');xlabel('L1');
title('L1 hopf\_branch');
subplot(3,3,[7,8]);
plot(a21_h2,hopf2testfuncs.genh(1,:),lw{:});
set(gca,'ylim',[-0.1,0],'xlim',a21lim);
grid on
xlabel('a_{21}');ylabel('L1');
title('First Lyapunov coefficient L1 along hopf2\_branch');
%% Periodic orbit continuation
% The convenience function |SetupPsol| creates the initial branch of
% periodic orbits. Its first arguments are |funcs|, the Hopf or
% steady-state branch along which a Hopf bifurcation was detected, and the
% index of the Hopf point along the branch where one wants to branch off
% point.
% Important additional parameters:
%
% * |'radius'|: amplitude of initial periodic orbit (default 1e-3),
% * |'contpar'|: index of continuation parameter (otherwise taken from input branch),
% * |'degree'|: degree of collocation polynomials (default 3)
% * |'intervals'|: number of collocation intervals (default 20)
% * |'submesh'|: (default |'linear'|) where polynomial is stored in each
% collocation interval. For high degrees one should use |'cheb'|.
% * |'collocation_parameters'|: (default |'legendre'|) at which points in
% each collocation interval will the differential equation be evaluated?
% Alternative: |'cheb'| (recommended for stability computation in neutral
% equations).
% * |'hopfcorrection'|: apply Newton iteration to find correct Hopf point
% before branching off. This is recommended if the point provided is
% not of kind |'hopf'|.
%
% All other parameters are passed on to the resulting output branch of
% periodic orbits (here |psol_branch|). The parameter field of the output
% branch is inherited from the input branch (before optional parameters are
% taken into account).
[psol_branch,suc]=SetupPsol(funcs,branch1,indhopf,'degree',4,'intervals',40,...
    'max_step',[indp.a21,0.1]);
figure(2); clf;ax2=gca;
[psol_branch,s,f,r]=br_contn(funcs,psol_branch,50,'plotaxis',ax2);
% look at the period along the branch:
figure(4); clf;
Plot2dBranch(psol_branch,'y',@(p)p.period,'funcs',funcs);
xlabel('a21');
ylabel('period');
%% Figure: family of periodic orbits
% Branch of periodic solutions emanating from a Hopf point. The
% branch turns at the far right.

%% Floquet multipliers of periodic orbits
% The function |GetStability| can be used for periodic orbits as well.
[nunst_per,dom,trivdef_per,psol_branch.point]=GetStability(psol_branch,...
    'funcs',funcs,'exclude_trivial',true);
indfold=find(nunst_per,1,'first');
fprintf('Loss of stability at point %d\n',indfold);
%% Continuation of folds of periodic orbits
% See <demo1_POfold.html> for extensive explanation. find inital
% approximate fold initialize problem functions and branch
[foldfuncs,pfold_branch]=SetupPOfold(funcs,psol_branch,indfold,...
    'contpar',[indp.a21,indp.taus],...
    'dir',indp.taus,'print_residual_info',1,'step',0.01,parbd{:});
% do not print residuals along branch
%branch5.method.point.print_residual_info=0;
% continue branch
figure(1);
pfold_branch=br_contn(foldfuncs,pfold_branch,100);
pfold_branch=br_rvers(pfold_branch);
pfold_branch=br_contn(foldfuncs,pfold_branch,60);
xlabel('a21');ylabel('tau_s');
title('Continuation of fold of periodic orbits');
%% Floquet multipliers of fold periodic orbits
[nunst_pf,dom,triv_defect,pfold_branch.point]=GetStability(...
    pfold_branch,'funcs',foldfuncs,'exclude_trivial',true);
fprintf('max number of unstable Floquet multipliers: %d\n',max(nunst_pf));
%% Defect of trivial Floquet multipliers
% plot the defect of the trivial Floquet multiplier as computed by
% GetStability
n_orbits=length(pfold_branch.point);
figure(5);clf;
a21_pfold=getp(pfold_branch,'a21');
taus_pfold=getp(pfold_branch,'taus');
plot(1:n_orbits,triv_defect,lw{:});
xlabel('point number');
ylabel('defect of trivial Floquet multiplier');
title('Comparison of Floquet multiplier computation with extended system');
%% Add fold of periodic orbits to Two-parameter bifurcation diagram
Plot2dBranch(pfold_branch,'funcs',foldfuncs,args{:});
%% save data:
save('demo1_simple_results.mat');
