%% Background information on defining systems and their solution
% The information below is not necessary for ordinary usage, but may help
% to find sources of errors (for example, if |p_correc| or |br_contn| do
% not succeed, or if one wants to use a general-purpose nonlinear solver
% routine (such as |fsolve|) or continuation package (such as |coco|).
%
% <html>
% $Id: demo1_definingsystems.m 367 2019-07-14 21:56:46Z jansieber $
% </html>
%% Set path and format
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...
    [base,'ddebiftool_utilities/']);
clear;
format compact
format short g
%#ok<*ASGLU,*NOPTS,*NASGU>
%% Example: Defining system for Hopf bifurcation
% During the demo |neuron| we encounter a Hopf bifurcation. Below we define
% the system and choose a rough estimate that can be used as initial
% guess for the Newton iteration.
parnames={'kappa','beta','a12','a21','tau1','tau2','taus'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
funcs=set_symfuncs(@sym_neuron,'sys_tau',@()[ip.tau1,ip.tau2,ip.taus]);
par0([ip.kappa,ip.beta,ip.a12,ip.a21,ip.tau1,ip.tau2,ip.taus])=...
     [    1/2,      -1,      1,  0.8,       0.2      0.2,     1.5];
hopf_est=p_tohopf(dde_stst_create('parameter',par0,'x',[0;0]),'funcs',funcs)
%% Right-hand side in the Newton iteration
% The function |p_correc| has a standard Newton iteration inside
% (|dde_nsolve|). This Newton iteration solves the defining system for each
% point type (according to the field |'kind'|) listed in the manual. The
% function defining the nonlinear system 0=f(u) is
%
%   function [res,J]=p_correc_rhs(funcs,method,point,free_par,...)
%
% For illustration, here is the right-hand side f for the above Hopf
% problem. The argument |hopf_est| serves as a template and indicates which
% problem is solved (according to |'kind'| the function |dde_kind_jac_res|
% is called) and free parameter |ip.a21|. Before arguments are passed on to
% |p_correc_rhs| they are preprocessed if indicated by the field
% |preprocess| in |method.point|.
hopf_mth=getfield(df_mthod('hopf'),'point');
hopf_mth.print_residual_info=1;
f0=@(x)p_correc_rhs(funcs,hopf_mth,hopf_est,ip.a21,'x',x)
%% Conversion between point structure and vector variable
% Furthermore, the point structure |hopf_est| is converted to a vector
% variable using |dde_x_from_point|, where the indices of free parameters
% are passed as a second argument. The routine |dde_ind_from_point| shows
% which component of the point structure is put where in the variable
% vector.
x0=dde_x_from_point(hopf_est,ip.a21)
indices=dde_ind_from_point(hopf_est,ip.a21)
disp('indices.x')
indices.x
disp('[indices.v.re, indices.v.im]')
[indices.v.re,indices.v.im]
%% Initial residual and Jacobian
% We observe that the system is not "square". It has one fewer equation
% than variables, since the last equation of the defining system in the
% manual is not actually enforced. 
[r0,J0]=f0(x0)
%% Newton iteration
% One may still pass this into a Newton solver, which will warn about
% non-square Jacobian. The general function
% |dde_point_from_x(x,template,free_par)| converts a vector back into
% point.
[xh0,suc]=dde_nsolve(f0,x0,hopf_mth)
hopf0_corrected=dde_point_from_x(xh0,hopf_est,ip.a21)
%% Appending extra conditions
% The function |p_correc_rhs| permits an optional argument |step_cond| that
% appends rows to the Jacobian (without enforcing the corresponding
% conditions). Thus, one might append the transpose of the nullvector of
% J0. Below the extended f is created and fed into the Newton iteration
% which does not complain anymore. The function |p_correc_rhs| has other
% optional arguments for passing on reference points (|pref|), or extra
% columns (|extracolumns|) to be appended to the Jacobian.
v0=dde_nullspace(J0); % same as null(J0) but also works for sparse J0
v0_pt=dde_point_from_x(v0,hopf_est,ip.a21) % convert vector back to point
f=@(x)p_correc_rhs(funcs,hopf_mth,hopf_est,ip.a21,'x',x,'step_cond',v0_pt)
[r,J]=f(x0);   % residual and Jacobian
size_J=size(J) % Jacobian is square now
[xh,suc]=dde_nsolve(f,x0,hopf_mth)
hopf_corrected=dde_point_from_x(xh,hopf_est,ip.a21)
%% Pre- and postprocessing inside |p_correc|
% The point methods contain fields |preprocess| and |postprocess| which are
% names of functions that are called inside |p_correc| before and after the
% Newton iteration to perform preprocessing (for example, creating
% additional problem specific rows for the Jacobian). The Hopf bifurcation
% continuation is an example for this.
hopf_mth
%% Setup for solution of nonlinear system
% The function |p_correc_setup| is a wrapper that calls |preprocess| if
% necessary and returns a function |f|, and initial guess |x0| and a struct
% |data| that can be passed on to the postprocessing. The following 4
% commands are the core of |p_correc| (except a recursive call for
% remeshing) in their general form (with additional optional arguments such
% as |'previous'| and |'step_cond'|). The function
%
% [r,J]=f(x0)
%
% should be "square", that is, it should contain as many variables as
% equations.
[f,x0,data]=p_correc_setup(funcs,hopf_est,ip.a21,hopf_mth)
[x,success]=dde_nsolve(f,x0,data.method)
data.point=dde_point_from_x(x,data.point,data.free_par)
hopf=dde_postprocess(data.point,data)
%% Other point types
% The routines |dde_x_from_point|, |dde_point_from_x|, |p_correc_rhs| work
% for all point types. A new point type |kind| needs to implement
% |dde_kind_create|, |dde_kind_jac_res| for the Newton iteration to work,
% and other functions such as |dde_kind_delays|, |dde_kind_dot|,
% |dde_kind_eig|, |dde_kind_delay_zero_cond| if needed (it can fall back
% onto predefined routines).