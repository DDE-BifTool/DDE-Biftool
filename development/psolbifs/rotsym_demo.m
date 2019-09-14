%% Demo of extension for rotational symmetry using Lang-Kobayashi equations
%
% The equation (modelling a laser subject to delayed coherent optical
% feedback) is given as
% 
% $$E'(t)=[1+i\alpha]n(t)E(t)+\eta\exp(i\phi)E(t-\tau), \qquad
%   n'(t)=\epsilon[p-n(t)-(2n(t)+1)\bar E(t)E(t)]$$
% 
% The main bifurcation parameters will be $\eta$ and $\phi$.
%
% *Warning* The functions for the extended systems determining relative
% equilibria (rotating waves) and relative periodic orbits (modulated
% waves), do not have support for user-provided system derivatves or
% state-dependent delays!
%
% <html>
% $Id: rotsym_demo.m 375 2019-09-14 14:26:50Z jansieber $
% </html>
%
%% Load DDE-Biftool and extension into Path
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...
    [base,'ddebiftool_extra_psol/'],...
    [base,'ddebiftool_utilities/'],...
    [base,'ddebiftool_extra_rotsym']);
load('../../demos/rotsym_demo/LKbifs.mat');
%% Problem definition using |set_rotfuncs|
% In addition to the user-defined functions |set_rotfuncs| needs the
% matrix |A| generating the rotation and (optional) the rotation as a
% function $\phi\mapsto \exp(A\phi)$. Then the system is assumed to have
% rotational symmetry $exp(A\phi)$ where |A| is anti-symmetric.
A=[0,-1,0; 1,0,0; 0,0,0];
expA=@(phi) [cos(phi), -sin(phi),0; sin(phi), cos(phi),0; 0,0,1];
%% Initial values of parameters and parameter indices
parnames={...
    'pump',...    % injection current
    'eta',...     % feedback strength
    'phi',...     % feedback phase
    'tau',...     %  delay
    'alpha',...   % alpha factor
    'epsilon',... % carrier relaxation time
    'omega'};     % rotation velocity
cs=[parnames;num2cell(1:length(parnames))];
ip=struct(cs{:});
par([ip.alpha,ip.pump,ip.epsilon,ip.eta,ip.phi,ip.tau])=...
    [    4,      0.1,    5e-3,     5e-3,   0,     100];
%% Right-hand side and call to |set_rotfuncs|
rfuncs=set_rotsymfuncs(@sym_LangKobayashi,...
    'rotation',A,'exp_rotation',expA,'sys_tau',@()ip.tau);
%% Fold of relative POs
% Folds of RPOs require a modified initialization routine since the
% condition |rot_cond| needs to be applied to the derivative, too. This
% requires again the introduction of an artificial parameter (automatically
% appended to parameter vector).  Theuser has to add rotation speed to the
% list of continuation parameters.
ind_mwf=find(diff(rw_phas_per_nunst)==1&real(dom_per(1:end-1))>0)+1;
[pfoldfuncs,mwfoldbr]=SetupMWFold(rfuncs,rw_phas_per,ind_mwf(1),...
    'contpar',[ip.phi,ip.eta,ip.omega],opt_inputs{:},...
    'print_residual_info',1,'dir',ip.eta);
%%
mwfoldbr=br_contn(pfoldfuncs,mwfoldbr,200,'plotaxis',ax);
mwfoldbr=br_rvers(mwfoldbr);
mwfoldbr=br_contn(pfoldfuncs,mwfoldbr,200,'plotaxis',ax);
%% Stability of orbits at fold of RPOs
% Here the stability includes the Floquet multiplier -1. Alternatively
% include  another 1 into the list in |'locate_trivial'|:
% |'locate_trivial',@(p)[1,1,1]|.
mwforbs=pfoldfuncs.get_comp(mwfoldbr.point,'solution');
[nunst_mwf,dom,triv_defect,mwforbs]=GetStability(mwforbs,...
    'exclude_trivial',true,'locate_trivial',@(p)[1,1],'funcs',rfuncs);
fprintf('max error of Floquet mult close to 1: %g\n',max(abs(dom-1)));
%%
save('LKbifs.mat');
%% re-plot all bifurcations
plot_2dbifs;
