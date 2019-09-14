%% test NDDE modifications - laser model by Hessel et al
%
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...
    [base,'ddebiftool_extra_psol/'],...
    [base,'ddebiftool_utilities/'],...
    [base,'ddebiftool_extra_rotsym/']);
clear;
format compact
format short g
lw={'linewidth',2};
%% Problem definition using |set_rotfuncs|
% In addition to the user-defined functions, |set_rotfuncs| needs the
% matrix |A| generating the rotation and (optional) the rotation as a
% function $\phi\mapsto \exp(A\phi)$. Then the system is assumed to have
% rotational symmetry $exp(A\phi)$ where |A| is anti-symmetric.
Crot=[0,-1;1,0];
Crotphi=@(phi)[cos(phi), -sin(phi); sin(phi), cos(phi)];
A=blkdiag(Crot,Crot);
expA=@(phi) blkdiag(Crotphi(phi),Crotphi(phi));
%% Set number of delays and create parameter names as strings
parnames={'a','b','h','tau','omega'};
cs=[parnames;num2cell(1:length(parnames))];
ip=struct(cs{:});
getp=@(br,name)arrayfun(@(p)p.parameter(ip.(name)),br.point);
%% Initial parameter and state values 
par0([ip.a,ip.b,ip.h,ip.tau,ip.omega])=...
     [   1, 0.0,  2,   5,    -1];
% par: h a b \Omega \Tau
%     [2 1 0  1      1]
x0=[1;0;0;0]; % initial state
%% Definition of problem
% rfuncs contains right-hand side and its derivative, rotation matrix, and
% left-hand mass matrix |M|. The mass matrix |M| is assumed to be constant.
% For systems with rotational symmetry the mass matrix must commute with
% the rotation matrix |A|: |AM=MA|. We also set parameter bounds and
% collect typical call arguments. Note that |rot_pointtype_list()| defines
% the known point types for problems with rotational symmetry. This is
% useful for stability computations, when detecting changes of stability
% (it helps ignoring trivial eigenvalues).
sfuncs=set_symfuncs(@sym_ddae_demo,'sys_tau',@()ip.tau,'p_vectorized',false);
rfuncs=set_rotfuncs('sys_rhs',sfuncs.sys_rhs,'sys_dirderi',sfuncs.sys_dirderi,...
    'rotation',A,'exp_rotation',expA,...
    'sys_tau',@()ip.tau,'x_vectorized',true,'lhs_matrix',diag([1,1,0,0]));
bd={'max_bound',[ip.b,2;ip.a,10],'min_bound',[ip.b,0.0; ip.a,-0.5]};
opt_inputs=[bd,{'extra_condition',1,'print_residual_info',1}];
stab_inputs={'exclude_trivial',true,'pointtype_list',@rot_pointtype_list};
%% Rotating wave solutions
% initialize branch for continuation in |b| (and |omega|)
rw_branch=SetupStst(rfuncs,'contpar',[ip.b,ip.omega],...
    'x',x0,'parameter',par0,opt_inputs{:},...
    'max_step',[ip.b,0.1]);
%% Continue roatating waves in |b|
figure(1);clf;ax1=gca;xlabel(ax1,'b');ylabel(ax1,'omega');
rw_branch=br_contn(rfuncs,rw_branch,70,'plotaxis',ax1);
rw_branch=br_rvers(rw_branch);
rw_branch=br_contn(rfuncs,rw_branch,10,'plotaxis',ax1);
bvals=getp(rw_branch,'b');
%% Stability computations of equilibria
% use the Breda discretization with Chebyshev polynomials. They increase
% the mesh adaptively until 20 eigenvalues nearest to 0 are acccurate. 
[stat_nunst,dom,defect,rw_branch.point]=GetStability(rw_branch, ...
    'funcs',rfuncs,stab_inputs{:},'nearest',0);
%% A simple animation of the spectrum
animate_spectrum(rw_branch,ip,'order','reverse','hold','on',...
    'pause',0.1,'figures',[3,4])
%% RW Hopf bifurcations
% This follows the standard procedure, except that |omega| needs to be
% specified as extra parameter. As expected the Hopf bifurcation point
% acquires other unstable eigenvalues at the line of essential instability
% |b=1|. Note that many unstable eigenvalues are invisible as they have
% high frequency.
ind_hopf=find(stat_nunst==0,1,'first');
[rw_hopf,suc]=SetupRWHopf(rfuncs,rw_branch,ind_hopf,'contpar',[ip.b,ip.a,ip.omega],...
    'dir',ip.b,'max_step',[0,0.1]);
figure(2);clf;ax2=gca;xlabel(ax2,'b');ylabel(ax2,'a');
rw_hopf=br_contn(rfuncs,rw_hopf,100,'plotaxis',ax2);
rw_hopf=br_rvers(rw_hopf);
rw_hopf=br_contn(rfuncs,rw_hopf,200,'plotaxis',ax2);
[hopf_nunst,dom,defect,rw_hopf.point]=GetStability(rw_hopf, ...
    'funcs',rfuncs,stab_inputs{:},'nearest',0);
disp('Hopf: b, a and number of unstable eigenvalues');
disp([getp(rw_hopf,'b');getp(rw_hopf,'a');hopf_nunst.']);
%% RW Fold bifurcations
% The RW fold setup creates a new functions structure as it solves an
% extended problem. Again |omega| needs to be specified as extra parameter.
% Since the fold has always |b>1| all folds are unstable in infinitely many
% directions. Our stability computation picks up the low frequency
% instabilities.
ind_fold=find(abs(diff(stat_nunst))==1,1,'first');
[foldfuncs,rw_fold,suc]=SetupRWFold(rfuncs,rw_branch,ind_fold,'contpar',[ip.b,ip.a,ip.omega],...
    'dir',ip.b);
figure(2);ax2=gca;
rw_fold=br_contn(foldfuncs,rw_fold,100,'plotaxis',ax2);
rw_fold=br_rvers(rw_fold);
rw_fold=br_contn(foldfuncs,rw_fold,100,'plotaxis',ax2);
[fold_nunst,dom,defect,rw_fold.point]=GetStability(rw_fold, ...
    'funcs',foldfuncs,stab_inputs{:},'nearest',0);
disp('Fold: b, a and number of unstable eigenvalues');
disp([getp(rw_fold,'b');getp(rw_fold,'a');fold_nunst.']);
%% Relative periodic orbits (modulated waves) branching off at Hopf
% Be aware that we use Chebyshev collocation (optional argument
% |'collocation_parameters'|), sparse matrix computation. We subdivide the
% period into 50 (adaptively adjusted) subintervals. On each subinterval
% the solution is approximately a polynomial, stored on the Chebyshev
% nodes. This non-equidistant mesh (argument |'submesh', 'cheb'|) permits
% arbitrary high degree, while standard mesh (default) is only safe for low
% degrees (<=8). The computation below could run with |'degree',30|, and
% fewer intervals. We stop slightly after the essential instability |b=1|.
degree=5;%
[mw_branch,suc]=SetupPsol(rfuncs, rw_branch,ind_hopf,...
    'intervals',50,'degree',degree,'submesh','cheb',...
    'collocation_parameters','cheb',...
    opt_inputs{:},'matrix','sparse','radius',0.1,'max_step',[0,0.1],...
    'max_bound',[ip.b,1.01])
figure(1);
mw_branch=br_contn(rfuncs,mw_branch,100,'plotaxis',ax1);
%% Stability of modulated waves
% The standard collocation (on Legendre points) appears to produce a large
% number of spurious unstable eigenvalues. The collocation using Chebyshev
% nodes 2 to degree+1 on each subinterval appears ok. This is not covered
% by any numerical analysis, though. The matrix size in this particular
% problem is only 804 (for intervals 50, degree 5, delay tau=5), such
% that|'eigmatrix', 'sparse'| may not be necessary. The call |eigs| may
% warn because it only obtains fewer eigenvalues than requested (default
% 20) to the desired accuracy.
%
% It is important to watch the eigenfunctions. If they show submesh
% oscillations, the corresponding Floquet multipliers are inaccurate.
[mw_nunst,dom,defect,mw_branch.point]=GetStability(mw_branch, ...
    'funcs',rfuncs,stab_inputs{:},'geteigenfuncs',true,'recompute',true,...
    'eigmatrix','sparse');
%% Animation of spectrum along branch
% We observe some torus bifurcations.
animate_spectrum(mw_branch,ip,'hold','off',...
    'pause',0.1,'figures',[3,4,5],'nunst',mw_nunst,'ip_essential','b')
%% Continuation of Torus bifurcation
% This follows standard procedure. We set the bound for b close to the
% essential instability boundary.
ind_tr=find(mw_nunst==2,1,'first');
[trfuncs,tr_branch,suc]=SetupMWTorusBifurcation(rfuncs,mw_branch,ind_tr,'contpar',[ip.b,ip.a,ip.omega],...
    'dir',ip.b,'print_residual_info',1,'max_bound',[ip.b,1.01;ip.a,10],...
    'max_step',[0,0.2])
figure(2);
tr_branch=br_contn(trfuncs,tr_branch,50,'plotaxis',ax2);
tr_branch=br_rvers(tr_branch);
tr_branch=br_contn(trfuncs,tr_branch,100,'plotaxis',ax2);
%% Stability of point at torus bifurcation
% The torus bifurcation converges onto the Hopf bifurcation close to the
% essential instability. Otherwise, it forms a stability boundary for most
% of its part. At some point it crosses a fold.
[tr_nunst,tr_dom,tr_def,tr_branch.point]=GetStability(tr_branch,...
    'funcs',trfuncs,stab_inputs{:},'geteigenfuncs',true,'eigmatrix','sparse');
%% Spectrum
% The animation of the spectrum confirms that the spectrum adn the
% bifurcation agree. We also see the fold crossing, when an eigenvalue
% crosses 1.
animate_spectrum(tr_branch,ip,'hold','off',...
    'pause',0.1,'figures',[3,4,5],'nunst',tr_nunst,'ip_essential','b')
%% Preparation of continuation of fold of modulated waves
% We extract the solution component where the stability changes in the
% Torus bifurcation. Then we use this solution to initialize the fold
% computation.
ind_pofold=find(tr_nunst==1,1,'first')-1;
mw2_branch=mw_branch;
foldini=trfuncs.get_comp(tr_branch,'solution');
foldini=foldini(ind_pofold);
mw2_branch.point=foldini;
[pofoldfuncs,pofold_branch,suc]=SetupMWFold(rfuncs,mw2_branch,1,'contpar',[ip.b,ip.a,ip.omega],...
    'dir',ip.b,'print_residual_info',1,'max_bound',[ip.b,1.01;ip.a,10],...
    'max_step',[0,0.2])
%% Continuation of folds of modulated waves
% These folds bend forwards and backwards but along part of the branch they
% form a stability boundary, as the stability computation below will
% reveal.
figure(2);
pofold_branch=br_contn(pofoldfuncs,pofold_branch,50,'plotaxis',ax2);
pofold_branch=br_rvers(pofold_branch);
pofold_branch=br_contn(pofoldfuncs,pofold_branch,200,'plotaxis',ax2);
%% Stability at the fold point
% with respect to other eigenvalues. This determines if the fold is a
% stability boundary.
[pofold_nunst,pofold_dom,pofold_def,pofold_branch.point]=...
    GetStability(pofold_branch,...
    'funcs',pofoldfuncs,stab_inputs{:},'geteigenfuncs',true,...'recompute',true,...
    'eigmatrix','sparse');
%%
animate_spectrum(pofold_branch,ip,'hold','off',...
    'pause',0.1,'figures',[3,4,5],'nunst',pofold_nunst,'ip_essential','b')
%% One parameter bifurcation diagram
% for |a=1|.
figure(1);clf;ax6=gca;hold(ax6,'on');
inp={'ax',ax6,stab_inputs{:},'max_nunst',4};
Plot2dBranch(rw_branch,'funcs',rfuncs,inp{:});
Plot2dBranch(mw_branch,'funcs',rfuncs,inp{:});
xlabel(ax6,'b');ylabel(ax6,'omega');
%% Two parameter bifurcation diagram
figure(2);clf;ax7=gca;hold(ax7,'on');
inp={'ax',ax7,stab_inputs{:},'max_nunst',2};
Plot2dBranch(rw_fold,'funcs',foldfuncs,inp{:});
Plot2dBranch(rw_hopf,inp{:});
Plot2dBranch(tr_branch,'funcs',trfuncs,inp{:});
Plot2dBranch(pofold_branch,'funcs',pofoldfuncs,inp{:});
set(ax7,'xlim',[0.6,1.2],'ylim',[0,10])
xlabel(ax7,'b');ylabel(ax7,'a');
