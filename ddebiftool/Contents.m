% DDEBIFTOOL
%
% Files
%   auto_eqd             - find cumulative discretization error on given mesh
%   auto_msh             - create new mesh on [0,1] equidistributing the error
%   biftool              - 
%   br_contn             - extend DDE-BIFTOOL branch
%   br_measr             - function [col,lengths]=br_measr(branch,measure)
%   br_plot              - plot branch
%   br_recmp             - recompute branch partially at selected points
%   br_refin             - add points to branch by clicking left mouse button,quit with middle or right mouse button
%   br_rvers             - function t_branch=br_rvers(branch)
%   br_splot             - 
%   br_stabl             - compute stability information along branch
%   correct_ini          - set and correct first (two) initial points along branch
%   dde_set_options      - parses varargin and assigns fields of structure options
%   delay_zero_cond      - residual & derivative wrt x, T and parameter for tau_d_nr(tz)=tau'_dnr(tz)=0
%   df_brnch             - set up empty branch with default values
%   df_derit             - approximate derivatives of delays wrt state and parameter with finite difference
%   df_deriv             - approximate derivatives of r.h.s wrt state and parameter with finite differences
%   df_measr             - define solution measure for plotting
%   df_mthod             - create methed structure with default parameters for continuation, solving and stability
%   fold_jac             - Jacobian and residual of nonlinear system for fold
%   get_lms_ellipse_new  - function [a,b]=get_lms_ellipse_new(order,delta,ellipse_nb)
%   get_pts_h_new        - creates cloud of points indicating area where eigenvalue can be
%   get_S_h_row          - function [S_h_row,nL,rcM]=get_S_h_row(AA,tau,hh,alpha,beta,interp_order)
%   hcli_eva             - function [px]=hcli_eva(profile1,t,x,m);
%   hcli_jac             - Jacobian and residual for connecting orbits
%   help_stst_stabil     - function [mu,nL]=help_stst_stabil(AA,tau,hh,alpha,beta,interp_order)
%   hopf_jac             - Jacobian for Hopf problem
%   mu_to_lambda         - lambda=mu_to_lambda(mu,h)
%   mult_app             - find Floquet multipliers & modes of periodic orbit
%   mult_crit            - find Floquet multiplier closest to unit circle and its mode
%   mult_plt             - function mult_plt(mu)
%   p_axpy               - function p=p_axpy(a,x,y)
%   p_correc             - correct point using Newton iteration
%   p_dot                - Compute dot product of two periodic solutions (psol structs)
%   p_measur             - function sc_measure=p_measur(p,measure)
%   p_norm               - norm of point for distance measuring
%   p_normlz             - function normalized_p=p_normlz(p)
%   p_pplot              - function p_pplot(point,component,colour)
%   p_remesh             - function rm_point=p_remesh(point,new_degree,new_mesh)
%   p_secant             - compute normalized secant
%   p_splot              - plot spectrum of point with stability information
%   p_stabil             - compute stability information for point
%   p_tau                - evaluate (state-dependent) delays along orbit
%   p_tofold             - convert stst point to fold point
%   p_tohcli             - convert point to connecting orbit
%   p_tohopf             - convert point to Hopf bifurcation point
%   p_topsol             - create starting guess for periodic solution derived from point
%   p_tostst             - convert point to stst (from fold, hopf, hcli or near branching)
%   p_tsgn               - find negative delays
%   poly_d2l             - values of second derivative of lagrange polynomials through 0:1/m:1 at c
%   poly_del             - values of derivative of lagrange polynomials through 0:1/m:1 at c
%   poly_dla             - values of derivative of lagrange polynomials through t at c
%   poly_elg             - return value(s) of Lagrange interpolation polynomial at c in [0,1]
%   poly_gau             - function c=poly_gau(m)
%   poly_lgr             - values of lagrange polynomials through t at c
%   poly_lob             - function c=poly_lob(m)
%   psol_eva             - function [px,it,cscal]=psol_eva(profile1,t,x,m);
%   psol_jac             - residual & Jacobian of collocation problem for periodic orbits
%   psol_msh             - create new mesh equidistributing the error
%   psol_sysvals         - call right-hand sides and derivatives and derivatives of delays in all points listed in xx
%   replace_branch_pars  - pass on optional arguments to branch structure
%   root_app             - approximate rightmost characteristic roots for equilibrium
%   root_cha             - Characteristic matrix and its derivatives
%   root_int             - first n rows of integration operator S(h) in R^(n x L)
%   root_nwt             - correct eigenvalues for equilibrium using Newton iteration
%   root_plt             - function root_plt(l0,l1,n1)
%   set_funcs            - fill in funcs structure with user-defined functions for use with DDE-Biftool functions
%   sparse_blkdiag       - sparse_blkdiag: create blockdiagonal sparse matrices
%   stst_jac             - residual and jacobian for equilibrium problem
%   stst_stabil          - compute spectrum of linearized r.h.s. in equilibrium
%   stst_stabil_nwt_corr - function stability=stst_stabil_nwt_corr(stability,AA,tau,method)
%   time_h               - recommend steplength for stability calculation (relative to max(tau))
%   time_lms             - function [alpha,beta]=time_lms(method,order)
%   time_nrd             - function gamma=time_nrd(epsi,r,s)
%   time_saf             - compute safety factor for given distance to imaginary axis 
%   VAopX                - vectorized matrix multiplication/linear system solving
