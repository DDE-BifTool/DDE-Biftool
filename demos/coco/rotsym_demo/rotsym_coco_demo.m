%% Demo of extension for rotational symmetry using coco (Lang-Kobayashi equations)
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
% waves), do not have support for user-provided system derivatves greater
% than 1 or state-dependent delays!
%
% <html>
% $Id$
% </html>
%
%% Load DDE-Biftool and extension into Path (make sure that coco paths are set too)
clear
base=[pwd(),'/../../../'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_extra_rotsym/'],...
    [base,'ddebiftool_utilities'],...
    [base,'ddebiftool_coco']);
format compact
format short g
%#ok<*NOPTS> 
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
    'omega'};     % rotation velocity (has to be last parameter!)
cs=[parnames;num2cell(1:length(parnames))];
ip=struct(cs{:});
par([ip.alpha,ip.pump,ip.epsilon,ip.eta,ip.phi,ip.tau])=...
    [    4,      0.1,    5e-3,     5e-3,   0,     100];
%% Right-hand side and call to |set_symrotfuncs|
% If the derivative of the right-hand side was symbolically generated, one
% may call |set_rotsymfuncs|, instead of |set_rotfuncs|. The delays and the
% rotation matrices and functions have to be specified separately.
rfuncs=set_rotsymfuncs(@sym_LangKobayashi,...
    'rotation',A,'exp_rotation',expA,'sys_tau',@()ip.tau);
%% Extended functions in rotating coordinates for rotating waves
% A rotating wave (relative equilibrium, RE) is a solution of the form
%
% $$ x(t)=\exp(A\omega t)x_0. $$
%
% The extended functions in |rfuncs| will always treat the user-provided
% system in rotating coordinates.
%
% Original system (by user):
%
% $$x'(t)=f(x(t),x(t-\tau_1),...,x(t-\tau_n))$$
%
% Rotating coordinates: $x(t)=\exp(A\omega t)y(t)$:
%
% $$y'(t)=-A \omega y(t)+f(y(t),\exp(-A\omega\tau_1)y(t-\tau_1),\ldots,
%    \exp(-A\omega\tau_n)y(t-\tau_n))$$
%
% The rotation speed is chosen such that the rotating wave
% $x(t)=\exp(A\omega t)x_0$ is turned into an equilibrium $y(t)=y_0$. This
% is achieved by solving for equilibria of the $y$ equation and adding a
% |sys_cond| (file |rot_cond.m|) to povide an equation determining
% $\omega$. For rotating waves this is
%
% $$ y_\mathrm{ref}^TAy=0 $$
%
% DDE-Biftool can now give the user's |sys_cond| access to reference
% points. This is controlled by the optional argument
% |'sys_cond_reference'| of |set_...funcs|. It is enabled for the functions
% created by |set_rot..funcs|.
%% Initial guess
% The extension for rotating and modulated waves is *not* able to cope with
% invariant equilibria. Thus, we generate a non-trivial rotating wave as
% our initial guess. For the laser the rotating waves correspond to
% stationary lasing (on state) $E(t)=E_0\exp(i\omega t)$, $n=n_0$.
[E0,n0,par]=LK_init(par,ip);
%% Bounds and monitors
% Coco has a mechanism to add monitors for continuation, which can then be
% bounded. Free parameters are always monitors.
bounds.eta=[0,1e-2];
bounds.phi=[-4,8];
rwplot_monitors={'power',@(p)norm(p.x(1:2))};
mwplot_monitors={...
    'powermax',@(p)max(sqrt(sum(p.profile(1:2,:).^2))),...
    'powermin',@(p)min(sqrt(sum(p.profile(1:2,:).^2))),...
    'power',@(p)sqrt(p_dot(p,p,'ind_comp',1:2))};
%% for stability computations 
% The |rotsym| extension has a different list of point types (different
% stability criteria). The stability monitoring constructor always needs
% this information.
add_stab=@(prob)dde_construct_stab(prob,'pointtypes',@rot_pointtype_list);
%% Relative equilibria varying phase |phi| of the delayed feedback
% The standard convenience function |SetupStst| works, but one must add the
% continuation parameter index |length(parameter)| to the index list of
% continuation parameters.
%
% *Warning* DDE-Biftool assumes that the rotation speed is the last
% parameter in the parameter vector!
%
% One may pass on the name omega as 'extpars' such that it always gets
% automatically added to the list of free parameters (and such that
% extra_condition is set to true). Similar to the other demos, for the
% initial |stst| branch the input |'info'| is the initial point. For some
% reason the nonlinear problem has difficulty converging in some points, so
% we adjust the corrector and continuation parameters slightly.
pt=dde_stst_create('x',[E0;n0],'parameter',par);
prob=coco_prob();
prob=dde_construct(prob,'stst','funcs',rfuncs,'info',pt,...
    'pnames',ip,'extpars','omega','monitors',rwplot_monitors);
prob=add_stab(prob);
prob=coco_set(prob,'cont','h_max',1e-1,'h_min',1e-5,'ItMX',800);
prob=coco_set(prob,'corr','ResTOL',1e-8);
% in the below call omega and power are not necessary, they ar efor output
% only
coco(prob, 'rw_phas', [], 1,{'phi','omega','power'},bounds.phi);
%% Extract stability and plot
bd=coco_bd_read('rw_phas');          % read run from file
rwbr=dde_read_branch('rw_phas');     % convert run output to branch as known from DDE-Biftool
nunst_eqs=coco_bd_col(bd,'USTAB')    % print number of unstable eigenvalues
figure(1);clf;ax1=gca;hold(ax1,'on');
lw={'linewidth',2};
txt={'Fontsize',16,'FontName','Courier','FontWeight','bold'};
set(ax1,lw{:},txt{:});
rwthm=dde_plot_theme('stst');
rwthm.bd.col2='power';
coco_plot_bd(rwthm,'rw_phas');
dde_plot_legend(ax1,'location','best');
%% Modulated waves (Relative periodic orbits, RPOs)
% Modulated waves are solutions of the form
%
% $$ x(t)=\exp(A\omega t)x_0(t)$$
%
% where $x0(t)=x0(t-T)$ for all $t$ and some period $T$. That is, $x(t)$ is
% quasi-periodic, but can be turned into a periodic solution in rotating
% coordinates. The transformation to rotating coordinates is the same as
% for relative equilibria, but the additional condition (|rot_cond|) is
%
% $$\int_0^1 y_\mathrm{ref}^TAy(t) \mathrm{d} t=0$$
%
% where $y_\mathrm{ref}(t)$ is a reference solution. Since DDE-Biftool does
% not give access to reference solutions in user-defined conditions,
% |rot_cond| returns residual 0 and Jacobian $y(t)^TA$.
%% RPOs branching off at 2nd Hopf of REs
% The initialization works with the standard routine. Again, the rotation
% speed needs to be added to the list of continuation parameters.
%% Branch off to periodic orbit
% We choose to branch off at Hopf point 3. For branching off at Hopf points
% toward periodic solutions, we avoid double- covering by setting
% |'bi_direct'| to |false|. We also prescribe the tangent (this is done by
% the constructor |dde_stst2spol|).
hdata=dde_read_solution('rw_phas','label','HB'); % get solution data at HB
prob=coco_prob();
prob=dde_construct(prob,'psol','data',hdata(3),'tangent',true,...
    'matrix','sparse','monitors',mwplot_monitors);
prob=coco_set(prob,'cont','NAdapt',5, 'bi_direct', false,...
    'h_max',10,'h_min',1e-4,'PtMX',200);
prob=coco_set(prob,'corr','SubItMX',1,'ItMX',5);
prob=add_stab(prob);
bdpsol=coco(prob, 'runpsol', [], 1,...
    {'phi','omega','powermax','powermin'},bounds.phi);
%% Plot
psbr=dde_read_branch('runpsol');
mwthm=dde_plot_theme('psol');
mwthm.bd.col2='powermax';
axes(ax1);coco_plot_bd(mwthm,'runpsol');
dde_plot_legend(ax1);
%% Plot of "phase portraits" of relative periodic orbits
figure(3);clf;ax3=gca;hold(ax3,'on');
Epow=@(x)sqrt(sum(x(1:2,:).^2,1));
for i=1:10:length(psbr.point)
    plot(ax3,Epow(psbr.point(i).profile),psbr.point(i).profile(3,:),'.-');
end
hold(ax3,'off');
grid(ax3,'on');
axis(ax3,'tight');
tdeco={'fontsize',14,'fontweight','bold'};
ldeco={'interpreter','LaTeX'};
set(ax3,tdeco{:});
xlabel(ax3,'$n$',ldeco{:});
ylabel(ax3,'$|E|$',ldeco{:});
%% Continuation of fold of relative equilibria
% The functions for the extended system are generated by |SetupRWFold|. The
% standard routine |SetupPOfold| has to be modified because the extended
% condition |rot_cond| has to be applied to the derivative, too. Also the
% rotation speed needs to be added to the list of continuation parameters.
% The extended system for fold continuation of REs has one additional
% artificial continuation parameter. Then label |rwfold| is necessary such
% that the constructor dde_stst2rwfold is chosen.
% 
fdata=dde_read_solution('rw_phas','label','SN');
nfolds=[1,3];
%%
figure(2);clf;ax2=gca;hold(ax2,'on');
for i=1:length(nfolds)
    run_name=sprintf('rwfold%d',i);
    prob=coco_prob();
    prob=dde_construct(prob,'rwfold','data',fdata(nfolds(i)),...
        'monitors',rwplot_monitors);
    prob=add_stab(prob);
    coco(prob, run_name, [], 1,{'phi','eta','omega','power'},...
        {bounds.phi,bounds.eta});
    rfthm=dde_plot_theme('fold');
    axes(ax2);coco_plot_bd(rfthm,run_name); %#ok<*LAXES>
    dde_plot_legend(ax2);
    drawnow
end
%% Replot from scratch in two parameters phi and eta
figure(2);clf;ax2=gca;hold(ax2,'on');
set(ax2,lw{:},txt{:});
clear bdfold foldbr
for i=length(nfolds):-1:1
    run_name=sprintf('rwfold%d',i);
    foldbr{i}=dde_read_branch(run_name);
    bdfold{i}=coco_bd_read(run_name);
    rfthm=dde_plot_theme('fold');
    axes(ax2);coco_plot_bd(rfthm,run_name);
    dde_plot_legend(ax2);
    drawnow
end
%% Continuation of Hopf bifurcations of relative equilibria
% For Hopf bifurcation continuation the standard routine |SetupHopf| works
% without modification (|SetupRWHopf| is a simple wrapper).
hdata=dde_read_solution('rw_phas','label','HB');
prob=coco_prob();
prob=dde_construct(prob,'hopf','data',hdata(2),...
    'monitors',rwplot_monitors);
prob=add_stab(prob);
coco(prob,'hopf', [], 1,{'phi','eta','omega','power'},...
    {bounds.phi,bounds.eta});
%%
hbr=dde_read_branch('hopf');
bdhopf=coco_bd_read('hopf');
rhthm=dde_plot_theme('hopf');
axes(ax2);coco_plot_bd(rhthm,'hopf');
dde_plot_legend(ax2);
%% Period doubling of relative POs
% The standard initialization works (wrapped to give it sensible name). The
% rotation speed needs to be added to the list of continuation parameters.
% The label is |'mwPD'|.
%%
pddata=dde_read_solution('runpsol','label','PD');
prob=coco_prob();
prob=dde_construct(prob,'mwPD','data',pddata(1),...
    'monitors',mwplot_monitors);
prob=coco_set(prob,'cont','NAdapt',1,'h_max',1,'h_min',1e-6);
prob=add_stab(prob);
pdfuncs=getfield(coco_get_func_data(prob,'dde','data'),'funcs');
coco(prob,'mwpd', [], 1,{'phi','eta','omega','powermax','powermin'},...
    {bounds.phi,bounds.eta});
%% Add period doubling to bifurcation diagram
bdpd=coco_bd_read('mwpd');
mpdthm=dde_plot_theme('PD');
axes(ax2);coco_plot_bd(mpdthm,'mwpd');
dde_plot_legend(ax2);
%% Check how far off trivial Floquet multipliers are and how many are unstable
% This is also stored in column |USTAB| of the bifurcation diagram.
pdbr=dde_read_branch('mwpd');
nunst_pd=coco_bd_col(bdpd,'USTAB');
pdorbs=pdfuncs.get_comp(pdbr,'solution_for_stability');
[~,~,triv_defect]=GetStability(pdorbs,'pointtype_list',@rot_pointtype_list,...
    'exclude_trivial',true);
fprintf('max error of Floquet mult close to -1: %g\n',max(abs(triv_defect)));
%% Fold of relative POs
% The label for initializing folds is 'mwPOfold'. Again (as in all previous
% runs) the parameter |omega| and the monitors do not have to be specified
% in the call to coco. The knowledge that |omega| needs to be free has been
% inherited from the previous runs.
%%
pfdata=dde_read_solution('runpsol','label','SN');
prob=coco_prob();
prob=dde_construct(prob,'mwPOfold','data',pfdata(1),...
    'monitors',mwplot_monitors);
prob=add_stab(prob);
prob=coco_set(prob,'cont','NAdapt',1,'h_max',1,'h_min',1e-6,'PtMX',200);
coco(prob,'mwpf', [], 1,{'phi','eta','omega','powermax','powermin'},...
    {bounds.phi,bounds.eta});
%% Add fold of periodic orbits to bifurcation diagram
mpfthm=dde_plot_theme('POfold');
axes(ax2);coco_plot_bd(mpfthm,'mwpf');
dde_plot_legend(ax2);
%% Check stability and error of trivial Floquet multipliers
bdpf=coco_bd_read('mwpf');
pfbr=dde_read_branch('mwpf');
pffuncs=getfield(coco_get_func_data(prob,'dde','data'),'funcs');
pforbs=pffuncs.get_comp(pfbr,'solution_for_stability');
nunst_pf=coco_bd_col(bdpf,'USTAB');
[~,dom,triv_defect]=GetStability(pforbs,'pointtype_list',@rot_pointtype_list,...
    'exclude_trivial',true);
fprintf('max error of Floquet mult close to 1: %g\n',max(abs(triv_defect)));
