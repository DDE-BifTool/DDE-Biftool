%% Normal form demo - bifurcation detection and normal form computation
% This file is self-contained and requires the extension
% |ddebiftool_extra_nmfm|.
%
% <html>
% $Id: nmfm_demo.m 109 2015-08-31 23:45:11Z jansieber $
% </html>
%% The system
% This system represents a non-dimensionalized model of two interacting
% layers of neurons.
% 
% $$ \dot{x}_1(t) = - x_1(t) - a g(bx_1(t-\tau_1)) + cg (dx_2(t-\tau_2)), $$
%
% $$ \dot{x}_2(t) = - x_2(t) - a g(bx_2(t-\tau_1)) + cg (dx_1(t-\tau_2)), $$
%
% where $g:\bf{R} \rightarrow \bf{R}$ is of the sigmoidal form
%
% $$ g(z) = \left[\tanh(z-1) + \tanh(1)\right]\cosh(1)^2. $$
% 
% The variables $x_1(t)$ and $x_2(t)$ represent the population-averaged
% neural activity at time $t$ in layers one and two, respectively.
% The parameter $a > 0$ is a measure of the strength of inhibitory feedback,
% while $c > 0$ measures the strength of the excitatory effect of one layer
% on the other. The parameters $b > 0$ and $d > 0$ are saturation rates and
% the delays $\tau_{1,2}$ represent time lags in the inhibitory feedback
% loop and excitatory inter-layer connection. Note that the system is
% symmetric with respect to interchanging the labels $1$ and $2$, so
% equilibria are necessarily of the form $(x_0,x_0)$.
%% Initialization
% First, we load the necessary DDE-Biftool extensions.
% We use existing files in the |funcs| structure. These files were
% generated gen_sys_nmfm_demo.,
clear variables;
close all;
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');
%% Set parameter names
% This could also be loaded from mat file |'FHN_parnames.mat'|.
parnames={'a','b','c','d','tau1','tau2'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Set funcs structure
funcs=set_symfuncs(@sym_nmfm_demo,'sys_tau',@() [ind.tau1, ind.tau2]);
%% Define the boundaries and step sizes of the active parameters.
n = 2;
astepsize = 0.0005;
cstepsize = 0.005;
amin = 0.0;
amax = 0.55;
cmin = 0.0;
cmax = 1.0;
%% Continuation and plot of steady state branch
% We set up a steady state branch using |SetupStst| and do a standard
% continuation. Then, we plot the branch.
branch1 = SetupStst(funcs,'x',[0;0],'parameter',[0.25, 2, 15/29, 1.2, 12.7, 20.2],...
    'contpar',1,'max_step',[1 astepsize],'max_bound',[1 amax],'min_bound',[1 amin],...
    'newheuristics_tests',0);
branch1.method.continuation.plot = 0;
[branch1,s,f,r] = br_contn(funcs,branch1,300);
branch1 = br_rvers(branch1);
[branch1,s,f,r] = br_contn(funcs,branch1,300);
[nunst_stst,dom,triv,branch1.point]=GetStability(branch1,'funcs',funcs);
figure;
[xm,ym] = df_measr(0,branch1,1);
br_plot(branch1,xm,ym);
xlabel('a');
ylabel('x_1');
axis([0 .5 -0.1 0.1])
%% Bifurcation detection on the steady state branch
% To detect bifurcations, we have to call the function |LocateSpecialPoints|.
%
% |function [branch,testfuncs,indices,biftype]=LocateSpecialPoints(funcs,branch,varargin)|
%
% Inputs:
%
% * |funcs|: problem functions
% * |branch|: branch of solutions (at the moment supported: hopf, fold, stst)
% * optional name-value pairs:
% * 'max_nferror' (default=|1e-3|) warn if difference between different
% approximation orders is larger than this (useful for finite-difference
% approximation only)
% * all other name-value pairs are passed on as fields for branch fields
%
% Outputs:
%
% * |branch|: updated branch with inserted special points
% * |testfuncs|: structure containing fields with branch specific test functions
% * |indices|: locations of special points in branch.point
% * |biftype|: cell array with names of bifurcations
%
% To adjust the behavior of this algorithm, you can change a
% number of parameters in |branch.method.bifurcation|:
%
% * |minimal_real_part|: only roots with a real part above this parameter
% will be taken into consideration when detecting bifurcations. The default
% is |-0.03|. Too many roots will cause overflow (because the hopf detector
% multiplies them all), so if an overflow is reported, try increasing this
% parameter.
% * |correction_tolerance|: to correct hopf bifurcations, the standard
% routine |p_correc| is used. The default value for this corrector is |1e-7|.
% * |secant_iterations|: codimension 2 bifurcations are corrected using a
% secant method. This is the maximum number of iterations, with default
% value 30.
% * |secant_tolerance|: this is the tolerance of the secant method, with
% default |1e-9|.
% * |radial_tolerance_factor|: after a bifurcation point has been corrected,
% it is checked whether the corrected point actually lies on the branch.
% If the bifurcation point originated from a change in the roots between
% points A and B, the new point is considered "on the branch" if its
% distance to the line AB is less than 25% of length(AB). This default
% percentage can be adjusted with this parameter.
%
% You will see some notifications about found hopf points not falling
% within the branch. This indicates a false positive: the correction
% algorithm produced a hopf point, but it ended up some distance away from
% the branch.
branch1.method.bifurcation.minimal_real_part = -0.03;
fprintf('Bifurcation detection\n');
[branch1,testfuncs,indices,biftype] = LocateSpecialPoints(funcs, branch1);
%% Extracting bifurcation locations
% As you can see, three Hopf points were detected. They will have been
% flagged, e.g. |branch1.point(132).flag = 'hopf'|. They will also have normal
% form information: |branch1.point(132).nmfm.L1 = -0.0133009607|. We will choose two
% hopf points and continue them. (These two were chosen because their branches
% happen to intersect, producing a Double Hopf bifurcation.)
%
% To this end, we need to know where the bifurcations were found. We do this by
% calling |br_getflags| on the steady state branch.
%
% |flagged_point_indices = br_getflags(branch)|
%
% Inputs:
%
% * |branch|: any branch
%
% Outputs:
%
% * |flagged_point_indices|: a matrix containing the indices of the
% bifurcation points in |branch|.
%
% The way the |FPI| (flagged point indices) matrix is structured is as follows. Every bifurcation has
% its own number (0 = |stst|, 1 = |hopf|, ...). You can convert between
% these representations by using |string = num2bif(number)| and
% |number = bif2num(string)|. Now, |FPI(bif2num('bifucation type'),:)| is a
% list of the indices at which the bifurcation of the designated type
% occurs. (For the full list, you can refer to any of these two files.)
%
% We first select the second Hopf point on the steady state branch.
FPI = br_getflags(branch1);
start_ind = FPI(bif2num('hopf'),2);
%% Continuation of the first Hopf branch
% We do a standard continuation, using the starting index obtained from the
% flagged point indices.
fprintf('----- Hopf branch 1 -----\n');
[branch2, suc] = SetupHopf(funcs, branch1, start_ind, 'contpar', [1 3],...
    'dir', 3, 'step', cstepsize,'plot',0);
branch2.parameter.min_bound=[5 0; 6 0; 3 cmin; 1 amin];
branch2.parameter.max_bound=[1 amax; 3 cmax];
branch2.parameter.max_step=[1 0.002; 3 cstepsize];
branch2.method.stability.minimal_time_step = 0.005; % default 0.01
[branch2,s,f,r]=br_contn(funcs,branch2,300);
branch2 = br_rvers(branch2);
[branch2,s,f,r]=br_contn(funcs,branch2,300);
%% Bifurcation detection on the first Hopf branch
% We repeat the steps for bifurcation detection, but now apply it to the
% hopf branch. This yields several higher-order bifurcations.
branch2.method.bifurcation.minimal_real_part = -0.03;
[branch2,testfuncs,indices,biftype] = LocateSpecialPoints(funcs, branch2);
%% The second Hopf branch
% We now repeat the entire exercise, but with the third hopf point as
% input.
fprintf('----- Hopf branch 2 -----\n');
FPI = br_getflags(branch1);
start_ind = FPI(bif2num('hopf'),3);
[branch3, suc] = SetupHopf(funcs, branch1, start_ind, 'contpar', [1 3], ...
    'dir', 3, 'step', cstepsize,'plot',0);
branch3.parameter.min_bound=[5 0; 6 0; 3 cmin; 1 amin];
branch3.parameter.max_bound=[1 amax; 3 cmax];
branch3.parameter.max_step=[1 0.002; 3 cstepsize];
branch3.method.stability.minimal_time_step = 0.005; % default 0.01
[branch3,s,f,r]=br_contn(funcs,branch3,300);
branch3 = br_rvers(branch3);
[branch3,s,f,r]=br_contn(funcs,branch3,300);
branch3.method.bifurcation.minimal_real_part = -0.03;
[branch3,testfuncs,indices,biftype] = LocateSpecialPoints(funcs, branch3);
%% Plotting the bifurcations on the branches
% We now plot the branches. On the hopf branches, we want to plot the
% bifurcation points with a special marker. This is done using the function
% |Plot2dBranch|.
figure(2);clf;ax2=gca;
hold(ax2,'on');
branch1.parameter.free=[1 3];
lg=Plot2dBranch(branch1,'ax',ax2);
lg=Plot2dBranch(branch2,'ax',ax2,'oldlegend',lg);
lg=Plot2dBranch(branch3,'ax',ax2,'oldlegend',lg);
xlabel('a');
ylabel('c');