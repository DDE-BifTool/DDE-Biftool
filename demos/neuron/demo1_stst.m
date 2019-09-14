%% Continuation and stability of steady states (equilibria)
%
% <html>
% $Id: demo1_stst.m 367 2019-07-14 21:56:46Z jansieber $
% </html>
%
% Once the user-defined functions are prepared, DDE-Biftool can compute and
% continue equilibria of the DDE, and compute their linearized stability,
% thus detecting local bifurcations. This demo requires <demo1_funcs.html> to
% have run beforehand.
%%
%#ok<*ASGLU,*NOPTS,*NASGU>
%
%% Initial guess for steady state
% It is clear that the neuron DDE has a steady state solution
% $(x_1^*,x_2^*)=(0,0)$ for all values of the parameters. We define a first
% steady state solution using the parameter values $\kappa=0.5$,
% $\beta=-1$, $a_{12}=1$, $a_{21}=2.34$, $\tau_1=\tau_2=0.2$ and
% $\tau_s=1.5$. Remember that we chose |par| $=[\kappa, \beta, a_{12},
% a_{21},\tau_1,\tau_2, \tau_s]$. It is recommended to use the provided
% function dde_stst_create as it fills in automatically all other fields
% and puts them into a default order. Arguments can be given in any order.
stst=dde_stst_create('parameter',[1/2, -1, 1, 2.34, 0.2, 0.2, 1.5],...
    'x',[0;0]);

%% Linear stability of initial equilibrium
% We get default point method parameters and correct the point, which,
% being already a correct solution, remains unchanged. Computing and
% plotting stability of the corrected point reveals it has one unstable
% real mode, see figure.

method=df_mthod('stst'); % this uses new Chebyshev approximation to find eigenvalues
method.stability.minimal_real_part=-1
[stst,success]=p_correc(funcs,stst,[],[],method.point)
% compute its stability:
stst.stability=p_stabil(funcs,stst,method.stability)
figure(1); clf;
p_splot(stst); % plot its stability:
%% Ask for roots with more negative real part
% In both figures, approximations $(\times)$ and corrections $(*)$ are
% nearly indistinguishable. The stability comptued by Chebyshev
% interpolation has a field err that indicates the residual of the
% approximation, when checked using the characteristic matrix.
method.stability.minimal_real_part=-4;
method.stability.max_number_of_eigenvalues=1000;
stst.stability=p_stabil(funcs,stst,method.stability); % recompute stability:
figure(2); clf;
p_splot(stst); % replot stability
%% Figures: Spectrum of equilibrium
% Approximated $(\times)$ and corrected $(*)$ roots of the characteristic
% equation of neuron system at its steady state solution
% $(x_1^*,x_2^*)=(0,0)$. Real parts computed up to $\Re(\lambda)\geq
% -\frac{1}{\tau}$ (top), $\Re(\lambda)\geq -2$ (bottom).

%% Initialize branch of trivial equilibria
% (See also <demo1_simple.html> how the steps below can be accomplished
% with a single call to SetupStst.) We will use this point as a first point
% to compute a branch of steady state solutions. First, we obtain an empty
% branch with free parameter $a_{21}$, limited by $a_{21}\in[0,5]$ and
% $\Delta a_{21}\leq 0.2$ between points.

% get an empty branch with ind_a21 as a free parameter:
branch1=df_brnch(funcs,ind_a21,'stst')
branch1.parameter
branch1.parameter.min_bound
% set bounds for continuation parameter
branch1.parameter.min_bound(end+1,:)=[ind_a21 0];
branch1.parameter.max_bound(end+1,:)=[ind_a21 5];
branch1.parameter.max_step(end+1,:)=[ind_a21 0.2];
% use stst as a first branch point:
branch1.point=stst;

%%  Extend and continue branch of trivial equilibria
% To obtain a second starting point we change  parameter value $a_{21}$
% slightly and correct again.Because we know how the branch of steady state
% solutions continued in $a_{21}$ looks like (it is constant at
% $(x_1^*,x_2^*)=(0,0)$) we disable plotting during continuation by setting
% the corresponding continuation method parameter to zero.

stst.parameter(ind_a21)=stst.parameter(ind_a21)+0.1;
[stst,success]=p_correc(funcs,stst,[],[],method.point)
% use as a second branch point:
branch1.point(2)=stst;
branch1.method.continuation.plot=0;

% continue in one direction:
[branch1,s,f,r]=br_contn(funcs,branch1,100)
% turn the branch around:
branch1=br_rvers(branch1);
% continue in the other direction:
[branch1,s,f,r]=br_contn(funcs,branch1,100)

%% Stability of branch of equilibria
% During continuation, sixteen points were successfully computed ($s=16$)
% before the right boundary $a_{21}=5$ was hit (signalled by a warning). No
% corrections failed ($f=0$) and no computed points were later rejected
% ($r=0$). Reversing the order of the branch points allows to continue to
% the left.
%
% After obtaining suitable measure structures we plot the real part of the
% approximated and corrected roots of the characteristic equation along the
% branch, (see figure). Notice the strange behaviour (coinciding of several
% complex pairs of roots) at $a_{21}=0$. At this parameter value one of the
% couplings between the neurons is broken. In fact, for $a_{21}=0$, the
% evolution of the second component is independent of the evolution of the
% first.
branch1.method.stability.minimal_real_part=-2;
branch1=br_stabl(funcs,branch1,0,0);

% obtain suitable scalar measures to plot stability along branch:
[xm,ym]=df_measr(1,branch1)
figure(3); clf;
br_plot(branch1,xm,ym,'b'); % plot stability along branch:
ym.subfield='l0';
br_plot(branch1,xm,ym,'c');
plot([0 5],[0 0],'-.');
axis([0 5 -2 1.5]);
xlabel('a21');ylabel('\Re\lambda');
% plot stability versus point number:
figure(4); clf;
br_plot(branch1,[],ym,'b');
br_plot(branch1,[],ym,'b.');
plot([0 30],[0 0],'-.');
xlabel('point number along branch');ylabel('\Re(\lambda)');
%% Figures: Stability of equilibria
%
% <html><a name=stststability></a>
% </html>
%
% Real parts of the approximated (top) and corrected (top,bottom) roots of
% the characteristic equation versus $a_{21}$ (top) respectively the point
% number along the branch (bottom).
%% Save and continue with Hopf bifurcations: <demo1_hopf.html>
save('demo1_stst_results.mat');
