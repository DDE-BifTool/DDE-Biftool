% DDEBIFTOOL_EXTRA_PSOL
%
% Files
%   app_dir_deriv         - 2nd order approx directional derivative in ith column of x0
%   extract_from_POfold   - extract components from pfold solution branch or point array
%   extract_from_tr       - extract components from extended solution branch or point array
%   generate_tau_ext      - generate array of additional delays needed for extended DDE in periodic fold continuation
%   POfoldInit            - crude initial guess for fold of periodic orbits
%   SetupPeriodDoubling   - SetupPeriodDoubling - Initialize continuation of period doubling bifurcation
%   SetupPOfold           - SetupPOfold - Initialize continuation of folds of periodic orbits
%   SetupTorusBifurcation - SetupTorusBifurcation - Initialize continuation of torus or period doubling bifurcations
%   sys_cond_POfold       - constraints used for extended DDE in periodic fold continuation
%   sys_cond_TorusBif     - additional conditions for torus and period doubling bifurcation
%   sys_deri_POfold       - partial derivatives of r.h.s of extended DDE for fold of periodic orbits
%   sys_deri_TorusBif     - partial derivatives of r.h.s of extended DDE for torus or period doubling bifurcation
%   sys_dtau_SD_PObif     - partial derivative of sys_tau for state-dependent delays
%   sys_rhs_POfold        - rhs of extended DDE for fold of periodic orbits
%   sys_rhs_SD_POfold     - rhs for fold of periodic orbits in SD-DDEs
%   sys_rhs_SD_TorusBif   - rhs for torus bifurcation of periodic orbits in SD-DDEs
%   sys_rhs_TorusBif      - r.h.s. of extended DDE for torus and period doubling bifurcation
%   sys_tau_SD_PObif      - delays of extended systems for fold, torus & period doubling
%   tauSD_ext_ind         - generate array of indices counting additional delays needed for extended SD-DDEs
%   TorusInit             - crude initial guess for stat of toruis bifurcation from Floquet mode
