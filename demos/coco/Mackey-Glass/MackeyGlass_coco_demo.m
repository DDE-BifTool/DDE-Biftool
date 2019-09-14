%% DDE-Biftool demo Mackey-Glass Equation using coco
%
% The Mackey-Glass equation is given by
% 
% $$x'(t)=\beta \frac{x(t-\tau)}{1+x(t-\tau)^n}-\gamma x(t)$$
% 
% Parameters are (in this order) |beta|, |n|, |tau| (|gamma| is not part of
% parameter vector).
%
% <html>
% $Id$
% </html>
%
%% load DDE-Biftool into path
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_extra_nmfm/'],...
    [base,'ddebiftool_utilities'],...
    [base,'ddebiftool_coco']);
format compact
format short g
%#ok<*ASGLU>
%% Initial parameters and state
parnames={'beta','n','tau','gamma'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
par0([ip.gamma,ip.beta,ip.n,ip.tau])=...
    [     1.0,      2,   10,    0];
x0=(par0(ip.beta)-1)^(1/par0(ip.n));
bounds.tau=[0,2];
bounds.beta=[0,5];
%% Set user-defined functions
%  Using |beta|,|n|,|tau| and |gamma| as parameters. Different ways to
%  define the right-hand side: |fvec| uses finite differences for
%  derivatives in a vectorized manner. |fsingle| does the same but does not
%  employ vectorization (should be slower), |fsymbolic| uses the right-hand
%  side and derivatives created via symbolic toolbox (see
%  <gen_sym_MackeyGlass.html>).
rhs=@(x,xtau,beta,n)beta*xtau./(1+xtau.^n)-gamma*x;
fvec=set_funcs(...
    'sys_rhs',@(xx,p)rhs(xx(1,1,:),xx(1,2,:),p(ip.beta),p(ip.n)),...
    'sys_tau',@()ip.tau,'x_vectorized',true);
fsingle=set_funcs(...
    'sys_rhs',@(xx,p)rhs(xx(1,1,:),xx(1,2,:),p(ip.beta),p(ip.n)),...
    'sys_tau',@()ip.tau);
fsymbolic=set_symfuncs(@sym_MackeyGlass,'sys_tau',@()ip.tau);
%% Choose which problem definition should be tested
funcs=fsymbolic;
%% Initialization of branch of non-trivial equilibria
pt_ini=dde_stst_create('x',x0,'parameter',par0);
prob=coco_prob();
prob=dde_construct(prob,'stst','funcs',funcs,'info',pt_ini,'pnames',ip);
prob=coco_set(prob,'cont','h_max',0.3);
prob=dde_construct_stab(prob);
bd=coco(prob, 'eqrun', [], 1,{'tau','xnorm'},bounds.tau);
%% Extract stability and plot
% |coco_bd_plot| is a provided routine for quick plotting. Otherwise, the
% standard |coco_bd_...| tools (such as |coco_bd_read|, |coco_bd_col|) can
% extract information from the bifurcation diagram. By default the
% DDE-Biftool point structure is stored in column |'point'|, its stability
% (if computed) is stored in column |'stability'|. |dde_read_branch| reads
% these and one solution to construct a branch structure as used by
% DDE-Biftool. Here, |eqbr| is the trivial solution branch in the format
% of DDE-Biftool.
bd=coco_bd_read('eqrun');
trivbr=dde_read_branch('eqrun');
nunst_eqs=coco_bd_col(bd,'USTAB')
figure(1);clf;ax1=gca;hold(ax1,'on');
ltx={'Interpreter','LaTeX'};
lw={'linewidth',2};
txt={'Fontsize',16,'FontName','Courier','FontWeight','bold'};
set(ax1,lw{:},txt{:});
coco_plot_bd(dde_plot_theme('stst'),'eqrun');
dde_plot_legend(ax1);
%% Continue Hopf bifurcation in two parameters
hdata=dde_read_solution('eqrun','label','HB');
prob=coco_prob();
prob=dde_construct(prob,'hopf','data',hdata);
prob=dde_construct_stab(prob);
bdhopf1=coco(prob, 'runhopf', [], 1,{'tau','beta'},{bounds.tau,bounds.beta});
%% Plot Hopf bifurcation in tau beta plane
bdhopf1=coco_bd_read('runhopf');
nunst_h1=coco_bd_col(bdhopf1,'USTAB');
figure(2);clf;ax2=gca;hold(ax2,'on');
set(ax2,lw{:},txt{:});
coco_plot_bd(dde_plot_theme('hopf'),'runhopf');
dde_plot_legend(ax2);
%% Branch off at Hopf bifurcation
prob=coco_prob();
prob=dde_construct(prob,'psol','data',hdata,'tangent',true,'matrix','sparse',...
    'intervals',20,'degree',4);
prob=coco_set(prob,'cont','NAdapt',10, 'bi_direct', false,'h_max',10,'h_min',1e-4,'PtMX',400);
prob=dde_construct_stab(prob);
bdpsol=coco(prob, 'runpsol', [], 1,{'tau','xnorm'},bounds.tau);
%% extract branch and plot
psolbr=dde_read_branch('runpsol');
axes(ax1);coco_plot_bd(dde_plot_theme('psol'),'runpsol');
%% Find period doubling bifurcations in two parameters
pddata=dde_read_solution('runpsol','label','PD');
prob=coco_prob();
prob=dde_construct(prob,'PD','data',pddata);
prob=dde_construct_stab(prob);
bdpd=coco(prob, 'runpd', [], 1,{'tau','beta'},{bounds.tau,bounds.beta});
%% Plot
axes(ax2);coco_plot_bd(dde_plot_theme('PD'),'runpd');
dde_plot_legend(ax2);
%% Branch off at period doubling 
% (Solutions at far end get inaccurate.)
pddata=dde_read_solution('runpsol','label','PD');
prob=coco_prob();
prob=dde_construct(prob,'psol','data',pddata,'tangent',true,'matrix','sparse',...
    'biftype','PD');
prob=coco_set(prob,'cont','NAdapt',10, 'bi_direct', false,'h_max',10);
prob=dde_construct_stab(prob);
bdp2=coco(prob, 'runp2', [], 1,{'tau','xnorm'},bounds.tau);
%%
p2br=dde_read_branch('runp2');
axes(ax1);coco_plot_bd(dde_plot_theme('psol'),'runp2');
dde_plot_legend(ax1);
%% Find 2nd period doubling bifurcations in two parameters
pd2data=dde_read_solution('runp2','label','PD');
prob=coco_prob();
prob=dde_construct(prob,'PD','data',pd2data);
prob=coco_set(prob,'cont','NAdapt',10,'h_max',10);
prob=dde_construct_stab(prob);
bdpd2=coco(prob, 'runpd2', [], 1,{'tau','beta'},{bounds.tau,bounds.beta});
%% Plot
axes(ax2);coco_plot_bd(dde_plot_theme('PD'),'runpd2');
dde_plot_legend(ax2);
%% Branch off at 3rd period doubling 
pd2data=dde_read_solution('runp2','label','PD');
prob=coco_prob();
prob=dde_construct(prob,'psol','data',pd2data,'biftype','PD','radius',1e-4);
prob=dde_construct_stab(prob);
prob=coco_set(prob,'cont','NAdapt',10, 'bi_direct', false,'h_max',10);
bdp4=coco(prob, 'runp4', [], 1,{'tau','xnorm'},bounds.tau);
%%
p4br=dde_read_branch('runp4');
axes(ax1);coco_plot_bd(dde_plot_theme('psol'),'runp4');
dde_plot_legend(ax1);
%% Find 3rd period doubling bifurcations in two parameters
pd4data=dde_read_solution('runp4','label','PD');
prob=coco_prob();
prob=dde_construct(prob,'PD','data',pd4data);
prob=dde_construct_stab(prob);
bdpd4=coco(prob, 'runpd4', [], 1,{'tau','beta'},{bounds.tau,bounds.beta});
%% Plot
axes(ax2);coco_plot_bd(dde_plot_theme('PD'),'runpd4');
dde_plot_legend(ax2);
