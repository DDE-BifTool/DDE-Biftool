%% DDE-BIFTOOL  state-dependent delays sd-demo
% 
% This demo is an illustrative example, showing how to perform bifurcation
% analysis for a system with state-dependent delays.
%
% The demo shows
%
% * which functions the user has to provide and how to put them into the
% structure |funcs|
% * continuation of equilibria and their linear stability,
% * detection and continuation of Hopf bifurcations,
% * branching off from Hopf bifurcation and continuation of periodic orbits
% 
% <html>
% $Id: sd_demo.m 20 2014-04-11 19:27:33Z jan.sieber $
% </html>
% 
%
%% 

%% Differential equations
% The differential equations for this example are
%
% $$\begin{array}{l}
% \frac{\mathrm{d}}{\mathrm{d} t}x_1(t)=\frac{1}{p_1+x_2(t)}\left(1-p_2x_1(t)x_1(t-\tau_3)
% x_3(t-\tau_3)+p_3x_1(t-\tau_1)x_2(t-\tau_2)\right),\\
% \frac{\mathrm{d}}{\mathrm{d} t}x_2(t)=\frac{p_4 x_1(t)}{p_1+x_2(t)}+
%         p_5\tanh(x_2(t-\tau_5))-1,\\
% \frac{\mathrm{d}}{\mathrm{d} t}x_3(t)=p_6(x_2(t)-x_3(t))-p_7(x_1(t-\tau_6)-x_2(t-\tau_4))e^{-p_8 \tau_5},\\
% \frac{\mathrm{d}}{\mathrm{d} t}x_4(t)=x_1(t-\tau_4)e^{-p_1 \tau_5} -0.1,\\
% \frac{\mathrm{d}}{\mathrm{d} t}x_5(t)=3(x_1(t-\tau_2)-x_5(t))-p_9,
% \end{array} $$
%
% where $\tau_1$ and $\tau_2$ are constant delays and
%
% $$
% \begin{array}{l}
% \tau_3=2+p_5\tau_1x_2(t)x_2(t-\tau_1),\\
% \tau_4=1-\frac{1}{1+x_1(t)x_2(t-\tau_2)},\\
% \tau_5=x_4(t),\\
% \tau_6=x_5(t).
% \end{array}
% $$
%
% This system has five components $(x_1,\ldots,x_5)$, six delays
% $(\tau_1,\ldots,\tau_6)$ and eleven parameters $(p_1,\ldots,p_{11})$,
% where $p_{10}=\tau_1$ and $p_{11}=\tau_2$. 
%
%%
clear;                           % clear variables
%close all;                       % close figures
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...  % add ddebiftool folders to path
    [base,'ddebiftool_extra_nmfm'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_utilities']);
%#ok<*ASGLU,*NOPTS,*NASGU>
%% Specification of state-dependent delays
% The function |sys_tau| returns the value of the delay. *This is in
% contrast to* the definition of |sys_tau| for constant delays, where the
% delay's _position_ in the parameter list is returned. For
% state-dependent delays, the header of the function |sys_tau|, defining the
% delays, has the format (see <sd_tau.html> for the demo example):
%
%   function tau=sys_tau(k,xx,par)
%
% it defines the |k| th delay $\tau_k$, which is allowed to depend on
% |xx(:,1:k)| (the states $x(t)$,..., $x(t-\tau_{k-1})$) and the parameter
% |par|. This definition of delays permits the user to create arbitrary
% levels of nesting. Note that, when |sys_tau| is called with first
% argument |k|, its second argument |xx| has column dimension |k|.
%
% *Note* The order of the delays corresponds to the order
% in which they appear in |xx| as passed to
% the functions |sys_rhs| and |sys_deri|.
%
% Optionally, the user is encouraged to provide a function that  supplies
% derivatives of all delays with respect to the state and parameters. Its functionality is similar to the function
% |sys_deri|. Its header has the format 
%
%   function dtau=sys_dtau(delay_nr,xx,par,nx,np)
%
% Its format is similar to |sys_deri|. The result |dtau| is a scalar,
% vector or matrix of partial derivatives of the delay with number
% |delay_nr| which depends on the type of derivative requested via
% |nx| and |np|. See <sd_dtau.html> for an example.

%% The number of delays
% For systems with state-dependent delays DDE-Biftool requires a separate
% function that provides the number of delays:
%
%   function n=sys_ntau()
%
% returns the number |n| of delays. Accordingly, the argument |xx| of
% the right-hand side |sys_rhs| will have the column dimension |n+1|. The
% corresponding field in the array |funcs| defining the user functions is
% named |'sys_ntau'|.

%% Assignment of user function fields
% The functions can have arbitrary names or can be anonymous. The above
% names are the field names in the structure containing th euser functions.
% Assign these fields automatically using |set_funcs|. Note that the flag
% |tp_del| is set to 1 for state-dependent delays. This is determined
% inside |set_funcs| using a try-catch enclosed call |sys_tau()|.
% Derivatives are optional and will be approximated by finite differences
% if not provided.
ufuncs=set_funcs('sys_rhs',@sd_rhs,'sys_tau',@sd_tau,...
    'sys_ntau',@()6,'sys_deri',@sd_deri,'sys_dtau',@sd_dtau,'x_vectorized',true)
nfuncs=set_funcs('sys_rhs',@sd_rhs,'sys_tau',@sd_tau,...
    'sys_ntau',@()6,'x_vectorized',true)
%% Automatic generation of right-hand side and derivatives
% Alternatively, one may generate the code for the right-hand side, the
% delays and their derivatives using the symbolic toolbox. See
% <gen_sym_sd_demo.html>. Then the function definition structure can be
% defined using te wrapper |set_symfuncs|:
sfuncs=set_symfuncs('sym_sd_demo')
%% Select the r.h.s. to be tested
funcs=sfuncs;
%% Define parameter names
ntau=6;
parnames=[strcat('p',num2cell('1':'9')),{'tau1','tau2'}];
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
%% Continuation and stability of Equilibria for state-dependent delays
% Once the user-defined functions are prepared, DDE-Biftool can compute and
% continue equilibria of the DDE, and compute their linearized stability,
% thus detecting local bifurcations.
%% Initial guess for steady state
% We define a steady state solution using the parameter values listed in
% |p0| and an initial guess given as named argument |'x'| to |SetupStst|.
% For the first steady-state branch the continuation parameter will be |p5|
% and it will be varied between |-1| and |1|. Note that a |0| in the first
% entry of the argument |'max_step'| refers to the predictor stepsize.
p0([ip.p1,ip.p2,ip.p3,ip.p4,ip.p5,ip.p6,ip.p7,ip.p8,ip.p9,ip.tau1,ip.tau2])=...
   [  4.5, 0.04, -1.4,    6,-0.45,-0.01,    3,  0.3,  0.1,      1,    0.2];
[stbr,suc]=SetupStst(funcs,'x',[1.4; 1.5; -25; 0.6; 1.4],'parameter',p0,...
    'contpar',ip.p5,'step',-0.1,...
    'min_bound',[ip.p5,-1],'max_bound',[ip.p5,1],'max_step',[0,0.1])
%% Continuation of equilibria
% With two starting points and suitable method parameters we continue the
% branch |stbr| (with plotting) versus parameter $p_5$, see figure below.
% During continuation, 12 points were successfully computed before the
% state-dependent delay function $\tau_3$ crossed zero (signalled by a
% warning). The computed point with $\tau_3<0$ was not accepted. Instead,
% the point corresponding to $\tau_3=0$ was computed, see figure below. We
% check the value of $\tau_3$ at the last point in the branch using
% function |p_tau|.
%
% In similar cases, it might happen that the computed value of a delay is 
% a very small negative value. Because stability cannot be computed when
% there are negative delays, small negative delay values are automatically
% neglected when their value is larger than the value defined in the field
% |method.stability.delay_accuracy|.
figure(1); clf;
[stbr,s,f,r]=br_contn(funcs,stbr,50) % continue with plotting the branch
hold on
plot(stbr.point(end).parameter(ip.p5),stbr.point(end).x(1),'o','linewidth',2);
xlabel('p5');ylabel('x(1)');
fprintf('delay %d=%g\n',3,p_tau(funcs,stbr.point(end),3));
%% Figure: branch of equilibria
% Output of |br_contn|: predictions and corrections
% after computation of a branch of steady state solutions versus parameter
% $p_5$. |o| - the last computed point in the branch (corresponding to
% $\tau_3=0$)
%% Stability of equilibria and eigenvalues of linearization
% We compute the stability along the branch and after obtaining suitable
% measure structures we plot the real part of the corrected roots of the
% characteristic equation along the branch versus the point numbers, see
% figure below. From this figure it is not clear which real parts
% correspond to real roots respectively complex pairs of roots. We check
% the first point with unstable eigenvalues
[nunst_stbr,~,~,stbr.point]=GetStability(stbr,'funcs',funcs);
[xm,ym]=df_measr(1,stbr);
figure(2); clf;
br_plot(stbr,[],ym,'b.-');   % stability along branch versus point number
%br_plot(branch1,[],ym,'b.');
plot([0 10],[0 0],'-.');      % axis
grid on
unst1=find(nunst_stbr>0,1,'first');
stbr.point(unst1).stability.l0 % stability of point 5
xlabel('point number');ylabel('\Re(\lambda)');
%% Figure: Eigenvalues of linearization along branch
% Real parts of the corrected roots of the characteristic equation along
% the branch.
%% Detection and continuation of Hopf bifurcations
% The eigenvalues of the linearized system along branches of equilibria
% indicate potential bifurcations. In this demo complex conjugate pairs of
% eigenvalues cross the imaginary axis, corresponding to Hopf bifurcations.
% The convenience function |LocateSpecialPoints| refines the branch to
% locate the Hopf point more accurately and determines criticality (the
% first Lyapunov coefficient |L1|)
[stbr_wbifs,tests,stbif_ind,stbif_type]=LocateSpecialPoints(funcs,stbr)
%% Continuation of Hopf points in two parameters
% The demo will proceed to continue two of these Hopf bifurcations in two
% system parameters $p_2$ and $p_9$. We select the Hopf point and turn it
% into a Hopf bifurcation point type. Then we correct the Hopf-like point
% using appropriate method parameters and one free parameter ($p_5$).
[hopfbr,suc]=SetupHopf(funcs,stbr_wbifs,stbif_ind,'contpar',[ip.p2,ip.p9],...
    'step',0.1,'dir',ip.p9,...
    'min_bound',[ip.p2,-1; ip.p9,-1],'max_bound',[ip.p2,0.053; ip.p9,10],...
    'max_step',[ip.p2,1;ip.p9,1])
figure(3); clf;
[hopfbr,s,f,r]=br_contn(funcs,hopfbr,30); % continue with plotting Hopf branch
xlabel('p2');ylabel('p9');
%% Continuation and stability of periodic orbits
% DDE-biftool computes one-parameter families of periodic orbits
% (automatically determining their period) by solving periodic
% boundary-value problems approximately with collocation schemes. A
% typical starting point for a family of periodic orbits is a Hopf
% bifurcation.
%% Constructing an initial small-amplitude orbit near a Hopf bifurcation
% We use the Hopf point in the |stbr| branch to construct a small amplitude
% (|1e-1|) periodic solution on an equidistant mesh of |15| intervals with
% piecewise polynomial degree |3| (these are optional arguments of
% |SetupPsol|. If we plan to change the continuation parameter, we have to
% explicitly provide it as the optional parameter |'contpar'|.
[psolbr,suc]=SetupPsol(funcs,stbr_wbifs,stbif_ind,'contpar',ip.tau1,...
    'radius',0.1,'intervals',15,'degree',3,...
    'max_bound',[ip.tau1,10],'min_bound',[ip.tau1,0],'max_step',[ip.tau1,1e-2])
figure(4); clf;hold on
[psolbr,s,f,r]=br_contn(funcs,psolbr,20); % compute periodic solutions branch
point=psolbr.point(end);
p_ampl=max(point.profile(1,:))-min(point.profile(1,:));
plot(point.parameter(ip.tau1),p_ampl,'o');
xlabel('tau1');ylabel('max(x1)-min(x1)');
%% Stopping criterion
% As in the case of computing |stbr|, we have a warning, |BR_CONTN
% warning: delay number_3 becomes negative.| indicating that the delay
% function $\tau_3(t)$ became negative at some point(s) on the period
% interval of the computed solution during continuation of the branch. The
% periodic solution with $\tau_3(t)$ negative is not accepted as the branch
% point. Instead, the following algorithm is executed. First, using the
% solution with $\tau_3(t)$ negative and a mesh refinement, a time point
% |tz| is computed at which $\tau_3(t)$ reaches its minimum. Then, a
% periodic solution is computed under the conditions,
% $$
% \tau_3(\mathtt{tz})=0,\qquad\mathrm{d}\tau_3(\mathtt{tz})/\mathrm{d}t=0.
% $$
% We compute and plot the delay $\tau_3(t)$ on the mesh of representation
% points at the last accepted point in the branch, see figure below.
tau3fun=@(t)p_tau(funcs,point,3,t);
figure(5); clf;
plot(point.mesh,tau3fun(point.mesh),'+-','linewidth',2);
[tmin,taumin]=fminbnd(tau3fun,0,1)
hold on;
plot(tmin,taumin,'o','linewidth',2);
grid on
xlabel('t/period');ylabel('tau3');
%% Figure: Periodic orbit with state-dependent delay equal to 0
% $\tau_3(t/T)$ at the last computed point. Crosses indicate representation
% points of the mesh used. The value of |taumin| shows that $\tau_3(t)$ has
% its minimal value at a point |tmin| between two representation points.
% The function |tau3fun| uses polynomial interpolation with the degree
% given in the field |point.degree|.
%% A second family of periodic orbits - correction of initial orbit
% Now we use the last Hopf point in the |psol2br| to compute a
% branch of periodic solutions as a function of the parameter $p_1$, see
% figure below.
[psol2br,suc]=SetupPsol(funcs,hopfbr,length(hopfbr.point),'contpar',ip.p1,...
    'max_bound',[ip.p1,10],'min_bound',[ip.p1,0],'max_step',[ip.p1,1e-2])
figure(6); clf;hold on
[psol2br,s,f,r]=br_contn(funcs,psol2br,40); % compute periodic solutions branch
psol=psol2br.point(end);
tau6fun=@(t)p_tau(funcs,psol,6,t);
figure(7); clf;
plot(psol.mesh,tau6fun(psol.mesh),'+-','linewidth',2);
[t2min,tau6min]=fminbnd(tau6fun,0,1)
hold on;
plot(t2min,tau6min,'o','linewidth',2);
grid on
xlabel('t/period');ylabel('tau6=x5');
%% Figure: delay $\tau_6=x_5$ at last point of branch
% $\tau_6(t/T)$ at the last computed point. Dots indicate representation
% points of the mesh used.
%
% The minimal value of the delay $\tau_6$ is a negative value. 
% The stability of the corresponding solution is computed if this value is
% larger than the one defined in |method.stability.delay_accuracy|.
% Periodic orbits also have a trivial Flqouet multiplier. The function
% |GetStability| automatically excludes the trivial Floquet multiplier from
% the stability count if requested.
[nunst_psol2br,dom2,triverr2,psol2br.point]=GetStability(psol2br,...
    'funcs',funcs,'exclude_trivial',true)
figure(8); clf;
p_splot(psol2br.point(end));
axis image;
xlabel('\Re\mu'); ylabel('\Im\mu');
%% Figure: Floquet multipliers of last periodic orbit in |psol2br|
%% Save results, end of demo |sd_demo|
save('sd_demo.mat');
