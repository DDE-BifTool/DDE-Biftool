%% Definition of user functions
%
% <html>
% $Id: demo1_funcs.m 367 2019-07-14 21:56:46Z jansieber $
% </html>
%
% (Please load |ddebiftool| into path first, see <demo1.html>.)
% To define a system, the user should provide Matlab functions defining the
% right-hand side (in the example called |neuron_sys_rhs|) and the indices
% for the delays in the parameter vector (called |neuron_tau| in our
% example).
%% Right-hand side
% A function defining the right-hand side $f$:
% 
%   function y=sys_rhs(xx,par)
% 
% This function has two arguments, |xx| $\in R^{n\times (m+1)}$, which contains 
% the state variable(s) at the present and in the past,
% |xx| $=[x(t), x(t-\tau_1), \ldots, x(t-\tau_m)]$,
% |par| $\in R^{1\times p}$ which contains the parameters, |par| $=\eta$.
%
% For the example, this is ($f$ is called |neuron_sys_rhs| in our example)
%
% This simple demo manually numbers the parameters in the array as
% |par(1)...par(7)|. See <demo1_simple.html> for a safer way to keep
% track of parameter ordering in the parameter array.
neuron_sys_rhs=@(xx,par)[...
    -par(1)*xx(1,1)+par(2)*tanh(xx(1,4))+par(3)*tanh(xx(2,3));....
    -par(1)*xx(2,1)+par(2)*tanh(xx(2,4))+par(4)*tanh(xx(1,2))];

%% Delays
% The delays $\tau_i$, $i=1\ldots,m$ are considered to be part of the
% parameters ($\tau_i=\eta_{j(i)}$, $i=1,\ldots,m$). This is natural since
% the stability of steady solutions and the position and stability of
% periodic solutions depend  on the values of the delays. Furthermore
% delays can occur both as a 'physical' parameter and as delay, as in
% $\dot{x}=\tau x(t-\tau)$. From these inputs the right-hand side $f$ is
% evaluated at time $t$. For equations with constant delays DDE-Biftool
% determines which parameters are delays by calling an argument-less
% function of the type
%
%   function d=sys_tau()
%
% In the example we order the parameters as |par| $=[\kappa,
% \beta, a_{12}, a_{21},\tau_1,\tau_2, \tau_s]$. Thus, (giving it the name
% |neuron_tau|):
neuron_tau=@()[5,6,7];
ind_a21=4;  % used later for continuation
ind_taus=7; % used later for continuation

%% Jacobians of user-provided functions
% Optionally (recommended) the user may also specify the partial
% derivatives of the user-defined functions with respect to states, delayed
% states and parameters. For constant delays only the derivatives of $f$
% are required. They should be provided as a function of the form
%
%   function J=sys_deri(xx,par,nx,np,v)
%
% providing the partial derivatives of first and second order (see file
% <neuron_sys_deri.html> for details). Alternatives for providing
% derivatives are discussed below.
%
%% Definition of structure |funcs|
% Similar to standard functions such as |ode45| DDE-Biftool's routines have
% an argument that defines the right-hand side. Since DDE-Biftool needs
% several user-defined functions (sys_rhs, sys_tau, optionally sys_deri,
% sys_cond etc) these functions are collected in a structure |funcs|. This
% structure |funcs| is best set up by calling the DDE-Biftool routine
% |set_funcs| with a sequence of name-value pairs. Each name-value pair
% corresponds to a field in the structure. Fields that are not listed as
% arguments of set_funcs get replaced by a default if possible.
%
% Possible argument names are:
% 
% * |'sys_rhs'| (default |sys_rhs| if file |sys_rhs.m| present in folder):
%    right-hand side |sys_tau|
% * |'sys_tau'| (default |@()[]|): function defining delays
% * |'sys_deri'| (default |[]|): function defining partial derivatives of
%      |sys_rhs|, if absent, finite differences will be used
% * |'sys_dirderi'| (default |[]|): function defining directional
%     derivatives of |sys_rhs|. This is an alternative to providing
%     |sys_deri|. |sys_dirderi| may be a cell array, where entry
%     sys_dirderi{1}(xx,par,dx,dpar) is the first derivative in direction
%     (dx,dpar). Orders 1 and 2 are used for standdard bifurcation
%     analysis, higher orders (up to 5, only with dp=0) for normal form
%     computations
% * |'sys_mfderi'| (default |[]|) higher-order mixed derivatives with
%      respect to the state variables (one may use |sys_dirderi| instead).
% * |'sys_ntau'| (default 0, only needed for state-dependent delays) number
%    of delays
% * |'sys_cond'| (default |@dummy_cond|) function providing extra conditions
% * |'sys_dtau'| (default |[]|, only needed for state-dependent delays):
%      function defining partial derivatives of |sys_tau|
% * |'sys_dirdtau'| (default |[]|, only needed for state-dependent delays):
%      function defining directional derivatives of |sys_tau|. 
% * |'x_vectorized'| (logical, default false) set to true if |sys_rhs|,
%    |sys_deri| (or |sys_dirderi|, (|sys_tau| and |sys_dtau| for SD-DDEs) accept an argument
%    |xx| with three dimensions. For periodic-orbit computations the
%    function will be called with all collocation points simultaneously if
%    |x_vectorized| is true.
%
% Other fields are |tp_del| (true if delays are state-dependent),
% |sys_deri_provided| (true if user has provided |sys_deri|) and
% |sys_dtau_provided| (true if user has provided |sys_dtau|).
%% Original setup of problem
% The setup below is the originally recommended way of defining a
% right-hand side, delays and its first and second derivatives.
ftraditional=set_funcs(...
    'sys_rhs',neuron_sys_rhs,...
    'sys_tau',neuron_tau,...
    'sys_deri',@neuron_sys_deri) %#ok<NOPTS>
%% Numerical finite differences
% Instead one my rely on numerical finite differences for the derivatives.
fnumeric=set_funcs(...
    'sys_rhs',neuron_sys_rhs,...
    'sys_tau',neuron_tau) %#ok<NOPTS>
%% Directional derivatives
% The function |neuron_sys_deri| is quite cumbersome. Instead, one may
% define directional derivatives:
dtanh=@(x)(1-tanh(x)^2);
ddtanh=@(x)(tanh(x)^2-1)*2*tanh(x);
neuron_sys_dirderi{1}=@(x,p,dx,dp)[...
    -dp(1)*x(1,1)-p(1)*dx(1,1)+...
     dp(2)*tanh(x(1,4))+p(2)*dtanh(x(1,4))*dx(1,4)+...
     dp(3)*tanh(x(2,3))+p(3)*dtanh(x(2,3))*dx(2,3);...
    -dp(1)*x(2,1)-p(1)*dx(2,1)+...
     dp(2)*tanh(x(2,4))+p(2)*dtanh(x(2,4))*dx(2,4)+...
     dp(4)*tanh(x(1,2))+p(4)*dtanh(x(1,2))*dx(1,2)];
neuron_sys_dirderi{2}=@(x,p,dx,dp)[...
    -2*dp(1)*dx(1,1)+...
     2*dp(2)*dtanh(x(1,4))*dx(1,4)+p(2)*ddtanh(x(1,4))*dx(1,4)^2+...
     2*dp(3)*dtanh(x(2,3))*dx(2,3)+p(3)*ddtanh(x(2,3))*dx(2,3)^2; ...
    -2*dp(1)*dx(2,1)+...
     2*dp(2)*dtanh(x(2,4))*dx(2,4)+p(2)*ddtanh(x(2,4))*dx(2,4)^2+...
     2*dp(4)*dtanh(x(1,2))*dx(1,2)+p(4)*ddtanh(x(1,2))*dx(1,2)^2];
fdirectional=set_funcs(...
    'sys_rhs',neuron_sys_rhs,...
    'sys_tau',neuron_tau,...
    'sys_dirderi',neuron_sys_dirderi) %#ok<NOPTS>
%% Derivatives generated by symbolic toolbox
% If the symbolic toolbox (7.0 and greater) is available, one may generate
% derivatives and right-hand side automatically (see script
% <gen_sym_demo1.html>). Gnu Octave's package symbolic may also work. The
% output of the symbolic codee generation needs to be wrapped. The function
% |set_symfuncs| provides this (it calls |set_funcs| internally).
fsymbolic=set_symfuncs(@sym_neuron,'sys_tau',neuron_tau);
%% Vectorization of right-hand side
% For speed-up during periodic orbit computations, one can vectorize the
% provided functions, which make them faster. One may then set the flag
% |'x_vectorized'| in the call to |set_funcs|.  Code generated with the
% symbolic toolbox is automatically vectorized by default. The vectorized
% versions would look as follows for |sys_rhs| (and then using numerical
% approximations for the derivatives):
neuron_sys_rhs_vec=@(xx,par)[...
    -par(1)*xx(1,1,:)+par(2)*tanh(xx(1,4,:))+par(3)*tanh(xx(2,3,:));....
    -par(1)*xx(2,1,:)+par(2)*tanh(xx(2,4,:))+par(4)*tanh(xx(1,2,:))];
fnumeric_vec=set_funcs(...
    'sys_rhs',neuron_sys_rhs_vec,...
    'sys_tau',neuron_tau,...
    'x_vectorized',true) %#ok<NOPTS>
%% Vectorization of right-hand side
% The directional derivatives are vectorized similarly:
dtanh=@(x)(1-tanh(x).^2);
ddtanh=@(x)2*(tanh(x).^2-1).*tanh(x);
neuron_sys_dirderi_vec{1}=@(x,p,dx,dp)[...
    -dp(1)*x(1,1,:)-p(1)*dx(1,1,:)+...
     dp(2)*tanh(x(1,4,:))+p(2)*dtanh(x(1,4,:)).*dx(1,4,:)+...
     dp(3)*tanh(x(2,3,:))+p(3)*dtanh(x(2,3,:)).*dx(2,3,:);...
    -dp(1)*x(2,1,:)-p(1)*dx(2,1,:)+...
     dp(2)*tanh(x(2,4,:))+p(2)*dtanh(x(2,4,:)).*dx(2,4,:)+...
     dp(4)*tanh(x(1,2,:))+p(4)*dtanh(x(1,2,:)).*dx(1,2,:)];
neuron_sys_dirderi_vec{2}=@(x,p,dx,dp)[...
    -2*dp(1)*dx(1,1,:)+...
     2*dp(2)*dtanh(x(1,4,:)).*dx(1,4,:)+p(2)*ddtanh(x(1,4,:)).*dx(1,4,:).^2+...
     2*dp(3)*dtanh(x(2,3,:)).*dx(2,3,:)+p(3)*ddtanh(x(2,3,:)).*dx(2,3,:).^2; ...
    -2*dp(1)*dx(2,1,:)+...
     2*dp(2)*dtanh(x(2,4,:)).*dx(2,4,:)+p(2)*ddtanh(x(2,4,:)).*dx(2,4,:).^2+...
     2*dp(4)*dtanh(x(1,2,:)).*dx(1,2,:)+p(4)*ddtanh(x(1,2,:)).*dx(1,2,:).^2];
fdirectional_vec=set_funcs(...
    'sys_rhs',neuron_sys_rhs_vec,...
    'sys_tau',neuron_tau,...
    'sys_dirderi',neuron_sys_dirderi_vec,...
    'x_vectorized',true) %#ok<NOPTS>
%% Select one of the ways to define the problem for the demo
funcs=fsymbolic;
%% Save and continue to continuation and stability of steady states <demo1_stst.html>
save('demo1_funcs_results.mat');
