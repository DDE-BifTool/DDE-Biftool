%% DDE-BIFTOOL Minimal demo for using coco
%
% <html>
% $Id: minimal_coco_demo.m 346 2019-05-13 05:41:50Z jansieber $
% </html>
%
% This demo illustrates how to use the coco wrappers for tracking local
% bifurcations of equilibria and periodic orbits in DDEs with contant
% delay. This demo requires the extensions |ddebiftool_extra_psol| and the
% auxiliary functions in |ddebiftool_utilities| and |dde_coco|. The example
% is the Duffing oscillator with delayed feedback discussed in the
% large-delay limit by Yanchuk & Perlikowski in (PRE79,0462211,2009):
%
% $$x''(t)+d*x'(t)+a*x(t)+x^3+b*[x(t)-x(t-\tau)]=0$$
%
% The parameters are $(\tau,a,b,d)$ (used in this order in the |parameter|
% vector).

%% Define path and system
%
% Make sure that coco is intialized by running its |startup| script. Then
% we load the folder of DDE-Biftool and its extensions into the Matlab
% path, and define the right-hand side and the delays. In this case the
% right-hand side has been created using the symbolic toolbox such that we
% create the structure containing the user-defined functions using
% |set_symfuncs|.
% 
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities',...
    '../../ddebiftool_coco');
format compact
format short g
%#ok<*ASGLU>
%% Set parameter names
% This could also be loaded from mat file |'minimal_demo_parnames.mat'|, as
% it needs to be consistent with the dedfinitions in
% |gen_sym_minimal_demo.m|.
parnames={'tau','a','b','d'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Definition of right-hand side and delay parameters
% The simplest way to define the right-hand side is by creating a
% right-hand-side function (see manual) of the format |y=rhs(xx,p)| with
% |[n,ntaup1,nvec]=size(xx)|, and a 1 x np parameter array, where |size(y)|
% is |n x nvec|. The dimension |n| is the system size, |xx(:,k,:)| is the
% value of |xx| at time |t| minus |(k-1)|th delay (|xx(:,1,:)| is the value
% of |xx| at the current time. Enabling the function to be called with a 3d
% array argument for |xx| and returning a 2d array y, permits simultaneous
% (vectorized) calls of the function in periodic-orbit continuations for
% speed-up. The function can be created in a separate m file or as an
% anonymous function handle.
%
rhs=@(xx,p)[...
     xx(2,1,:);...
    -p(ind.d)*xx(2,1,:)-p(ind.a)*xx(1,1,:)-xx(1,1,:).^3-p(ind.b)*(xx(1,1,:)-xx(1,2,:))];
%% Call to set_funcs(...)
% For DDEs with parametric delays, the r.h.s. is passed on as named input
% |'sys_rhs'|, and the function indicating which parameters are delays is
% passed on as |'sys_tau'|. Set the optional input |'x_vectorized'| in
% |set_funcs(...)| to |true| to indicate vectorization of the r.h.s.
% function.
%funcs=set_funcs('sys_rhs',rhs,'sys_tau',@()ind.tau,'x_vectorized',true); %#ok<NASGU>
%% Definition of r.h.s. and derivatives using symbolic toolbox
% The function sym_minimal_demo has been created with the symbolc toolbox
% in <gen_sym_minimal_demo.html>. The wrapper |set_symfuncs| converts this
% function into a structure as needed for further DDE-Biftool computations.
% (|x_vectorized| and |p_vectorized| are on by default in this case.)
% Uncomment the statement below to use the simple definition of |rhs|
% above.
funcs=set_symfuncs(@sym_minimal_demo,'sys_tau',@()ind.tau);
%% Set bifurcation parameter range
bounds.tau=[0,20];
bounds.b=[0,0.61];
%% Definition of initial (trivial) equilibrium
% The system has only the trivial equilibrium for positive values of the
% parameter $a$. We create an example stst point (with dde_stst_create).
%#ok<*NOPTS>
par([ind.tau,ind.a,ind.b,ind.d])=...
    [    0,    0.5,  0.6, 0.05];
pt=dde_stst_create('x',[0;0],'parameter',par);
%% Trivial equilibrium branch
% The continuation below constructs the family of trivial equilibria.
% |dde_construct| contructs the branch by adding the zero problem to prob,
% using initial values from input |'info'| (also initial tangent, if
% requested). Names for parameters are passed on in |'pnames'|. Internally,
% the constructor also appends a monitoring function xnorm, which takes the
% L2 norm of the dynamic variables (field |'x'| or field |'profile'|).
%
% The contructor dde_construct_stab adds a monitoring function for
% stability changes. It also (by default) outputs the field stability into
% the bifurcation diagram.
prob=coco_prob();
prob=dde_construct(prob,'stst','funcs',funcs,'info',pt,'pnames',ind);
prob=dde_construct_stab(prob);
bd=coco(prob, 'runtriveqs', [], 1,{'tau','xnorm'},bounds.tau);
%% Extract stability and plot
% |coco_bd_plot| is a provided routine for quick plotting. Otherwise, the
% standard |coco_bd_...| tools (such as |coco_bd_read|, |coco_bd_col|) can
% extract information from the bifurcation diagram. By default the
% DDE-Biftool point structure is stored in column |'point'|, its stability
% (if computed) is stored in column |'stability'|. |dde_read_branch| reads
% these and one solution to construct a branch structure as used by
% DDE-Biftool. Here, |trivbr| is the trivial solution branch in the format
% of DDE-Biftool.
bd=coco_bd_read('runtriveqs');
trivbr=dde_read_branch('runtriveqs');
nunst_eqs=coco_bd_col(bd,'USTAB')
figure(1);clf;ax1=gca;hold(ax1,'on');
ltx={'Interpreter','LaTeX'};
lw={'linewidth',2};
txt={'Fontsize',16,'FontName','Courier','FontWeight','bold'};
set(ax1,lw{:},txt{:});
coco_plot_bd(dde_plot_theme('stst'),'runtriveqs');
%% Find labels of Hopf bifurcations
% The stability monitoring routine outputs labels such as |SN| and |HB|.
% Their label numbers can be found using standard coco tools.
ind_hopf=coco_bd_labs(bd,'HB');
%% Continue 2nd Hopf bifurcation
% |dde_read_solution| is a wrapper around |coco_read_solution|. Its output
% is an array of structs. Each of tem can be used to construct a hopf
% bifurcation problem in this case.
hdata=dde_read_solution('runtriveqs','label','HB');
prob=coco_prob();
prob=dde_construct(prob,'hopf','data',hdata(2));
prob=dde_construct_stab(prob);
bdhopf1=coco(prob, 'runhopf1', [], 1,{'tau','b'},{bounds.tau,bounds.b});
%% Plot
% A Hopf continuation detects double-Hopf points and Fold-Hopf points as HB
% and as HB.
bdhopf1=coco_bd_read('runhopf1');
nunst_h1=coco_bd_col(bdhopf1,'USTAB');
figure(2);clf;ax2=gca;hold(ax2,'on');
set(ax2,lw{:},txt{:});
coco_plot_bd(dde_plot_theme('hopf'),'runhopf1');
%% Double Hopf points
% We can branch off at second double Hopf bifurcation to the other Hopf
% bifurcation. The constructor for the emanating torus bifurcations is tbd
% for coco as tis requires normal form monitoring. (|C1branch_from_C2point| is te DDE-Biftool routine.)
hhdata=dde_read_solution('runhopf1','label','HB');
prob=coco_prob();
prob=dde_construct(prob,'hopf','data',hhdata(2));
prob=dde_construct_stab(prob);
bdhopf2=coco(prob, 'runhopf2', [], 1,{'tau','b'},{bounds.tau,bounds.b});
%% Plot
bdhopf2=coco_bd_read('runhopf2');
nunst_h2=coco_bd_col(bdhopf2,'USTAB');
axes(ax2);coco_plot_bd(dde_plot_theme('hopf'),'runhopf2');
%% Branch off to periodic orbit
prob=coco_prob();
prob=dde_construct(prob,'psol','data',hdata(2),'tangent',true,'matrix','sparse');
prob=coco_set(prob,'cont','NAdapt',20, 'bi_direct', false,'h_max',10,'h_min',1e-4,'PtMX',400);
prob=dde_construct_stab(prob);
bdpsol=coco(prob, 'runpsol', [], 1,{'tau','xnorm'},bounds.tau);
%% Plot
psbr=dde_read_branch('runpsol');
axes(ax1);coco_plot_bd(dde_plot_theme('psol'),'runpsol');
%% Continue Torus bifurcations
trdata=dde_read_solution('runpsol','label','TR');
clear bdtr
figure(3);clf;ax3=gca;hold(ax3,'on');
for i=length(trdata):-1:1
    rname=sprintf('runtr%d',i);
    prob=coco_prob();
    prob=dde_construct(prob,'torus','data',trdata(i),'matrix','sparse');
    prob=coco_set(prob,'cont','NAdapt',10,'h_max',10,'h_min',1e-4,'PtMX',400);
    prob=dde_construct_stab(prob);
    bdtr{i}=coco(prob, rname, [], 1,{'tau','b'},{bounds.tau,bounds.b});
    axes(ax3);
    coco_plot_bd(dde_plot_theme('torus'),rname);
    drawnow
end
%% Plot torusbifs
axes(ax2);
for i=length(trdata):-1:1
    rname=sprintf('runtr%d',i);
    bdtr{i}=coco_bd_read(rname);
    coco_plot_bd(dde_plot_theme('torus'),rname);
end
%% Continue POfold
pfdata=dde_read_solution('runpsol','label','SN');
prob=coco_prob();
prob=dde_construct(prob,'POfold','data',pfdata,'matrix','sparse');
prob=coco_set(prob,'cont','NAdapt',10,'h_max',10,'h_min',1e-4,'PtMX',400);
prob=dde_construct_stab(prob);
bdpf=coco(prob, 'runpf', [], 1,{'tau','b'},{bounds.tau,bounds.b});
%% Plot
bdpf=coco_bd_read('runpf');
axes(ax2);
coco_plot_bd(dde_plot_theme('POfold'),'runpf');
set(ax2,lw{:},txt{:});
set(ax1,lw{:},txt{:});
