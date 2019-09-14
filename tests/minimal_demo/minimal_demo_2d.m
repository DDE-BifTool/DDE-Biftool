%% Minimal demo - Equilibria, Hopf bifurcations, periodic orbits
%
% $Id: minimal_demo_2d.m 212 2017-07-09 22:31:00Z jansieber $
%
%#ok<*NOPTS>
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities',...
    '../tools');
format compact
format short g
%% Set number of delays and parameter names
ntau=1;
parnames={'tau','a','b','d'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Choose whether to test numerical differentiation
% Alternatively, one may set |'sys_mfderi',@minimal_demo_mfderi| instead of
% |'sys_dirderi',@minimal_demo_dirderi|. The function |minimal_demo_dirderi| was generated
% automatically using the maple script |maple_dirderi.mpl|, loading maple
% functions from |ddebiftool_deri.mpl|.
%
% The function sym_minimal_demo has been created with the symbolc toolbox
% in <gen_sym_minimal_demo.html>.
%
rhsxvec=@(x,p)[x(2,1,:);...
    -p(ind.d)*x(2,1,:)-p(ind.a)*x(1,1,:)-x(1,1,:).^3-p(ind.b)*(x(1,1,:)-x(1,2,:))];
fnumxvec=set_funcs('sys_rhs',rhsxvec, 'sys_tau',@()ind.tau,'x_vectorized',true);
drhsxvec={@(x,p,dx,dp)[dx(2,1,:);...
    -dp(ind.d)*x(2,1,:)-p(ind.d)*dx(2,1,:)-...
     dp(ind.a)*x(1,1,:)-p(ind.a)*dx(1,1,:)-...
     3*x(1,1,:).^2.*dx(1,1,:)-...
     dp(ind.b)*(x(1,1,:)-x(1,2,:))-p(ind.b)*(dx(1,1,:)-dx(1,2,:))]};
fdnumxvec=set_funcs('sys_rhs',rhsxvec, 'sys_tau',@()ind.tau,...
    'sys_dirderi',drhsxvec,'x_vectorized',true);
rhsxpvec=@(x,p)[x(2,1,:);...
    -p(1,ind.d,:).*x(2,1,:)-p(1,ind.a,:)*x(1,1,:)-x(1,1,:).^3-...
    p(1,ind.b,:).*(x(1,1,:)-x(1,2,:))];
fnumxpvec=set_funcs('sys_rhs',rhsxpvec, 'sys_tau',@()ind.tau,...
    'x_vectorized',true,'p_vectorized',true);
drhsxpvec={@(x,p,dx,dp)[dx(2,1,:);...
    -dp(1,ind.d,:).*x(2,1,:)-p(1,ind.d,:).*dx(2,1,:)-...
     dp(1,ind.a,:).*x(1,1,:)-p(1,ind.a,:).*dx(1,1,:)-...
     3*x(1,1,:).^2.*dx(1,1,:)-...
     dp(1,ind.b,:).*(x(1,1,:)-x(1,2,:))-p(1,ind.b,:).*(dx(1,1,:)-dx(1,2,:))]};
fdnumxpvec=set_funcs('sys_rhs',rhsxpvec, 'sys_tau',@()ind.tau,...
    'sys_dirderi',drhsxpvec,'x_vectorized',true,'p_vectorized',true);
rhs_simple=@(x,p)[x(2,1);...
    -p(ind.d)*x(2,1)-p(ind.a)*x(1,1)-x(1,1)^3-p(ind.b)*(x(1,1)-x(1,2))];
fsimple=set_funcs('sys_rhs',rhs_simple, 'sys_tau',@()ind.tau);
fsymbolic=set_symfuncs(@sym_minimal_demo,'sys_tau',@()ind.tau);
%% Choose one of |fnum|, |fsimple| or |fsymbolic|
fchoice={fsymbolic,fsimple,fnumxvec,fdnumxvec,fnumxpvec,fdnumxpvec};
%% set which r.h.s to be used
%
funcs=fchoice{indfuncs};
%% load reference solution
s=load('minimal_demo_2dselection.mat');
%% Trivial equilibria
s_e=@(s1,s2,tol)struct_error(rmfield(s1,'nvec'),rmfield(s2,'nvec'))<tol;
[triv_eqs,suc]=SetupStst(funcs,'x',[0;0],'parameter',[0.2,0.5,0.6,0.05],...
    'contpar',ind.tau,'max_step',[ind.tau,0.3],'max_bound',[ind.tau,20],...
    'newheuristics_tests',0,'minimal_real_part',-1);
assert(suc>0)
%% Stability of trivial equilibria
disp('Trivial equilibria');
figure(1);clf;ax1=gca;
[triv_eqs,succ,fail,rjct]=br_contn(funcs,triv_eqs,60,'plotaxis',ax1);
assert(succ>0&&rjct==0&&fail<succ/4)
nunst_eqs=GetStability(triv_eqs,'funcs',funcs);
%% Find Hopf bifurcations and their criticality
[triv_eqs,~,ind_hopf,stst_bifs]=LocateSpecialPoints(funcs,triv_eqs) 
%assert(s_e(triv_eqs.point(ind_hopf),s.triv_eqs.point(s.ind_hopf),1e-6));
%% Initialize 2nd Hopf bifurcation
hopfopts={'max_step',[ind.tau,0.2;ind.b,0.01],'max_bound',[ind.b,0.6]};
[hopf,suc]=SetupHopf(funcs,triv_eqs,ind_hopf(2),'contpar',[ind.b,ind.tau],'dir',ind.b,...
    hopfopts{:},'step',-1e-3);
assert(suc>0)
%% continuation of Hopf
figure(2);clf;ax2=gca;
[hopf,succ,fail,rjct]=br_contn(funcs,hopf,200,'plotaxis',ax2);
assert(succ>0&&rjct==0&&fail<succ/4)
%% Locate codim-2 bifurcations
[hopf_wbifs,htestfuncs,ind_hc2,hc2bifs]=LocateSpecialPoints(funcs,hopf);
%% test genh
genh=hopf_wbifs.point(ind_hc2(strcmp(hc2bifs,'genh')));
hohos=hopf_wbifs.point(ind_hc2(strcmp(hc2bifs,'hoho')));
%assert(s_e(genh,s.genh,1e-4));
%% test hoho
%assert(s_e(hohos,s.hohos,1e-6));
%% Branch off at genh to POfold
% We also find the first hopf bifurction of the branch |triv_eqs| and
% continue that, too.
[LPCfuncs,LPCbranch,suc]=C1branch_from_C2point(funcs,genh,hopf_wbifs.parameter.free,...
    'codim1','POfold','codim2','genh','print_residual_info',1,hopfopts{:});
assert(all(suc(:)>0))
%% continuation of LPC
[LPCbranch,succ,fail,rjct]=br_contn(LPCfuncs,LPCbranch,10,'plotaxis',ax2);
assert(succ>0&&rjct==0&&fail<succ/4)
%% Branch off to other Hopfs at Hopf-Hopf interactions
nhsteps=5;
hbr=hopf;
indh=1;
for i=length(hohos):-1:1
    hh=hohos(i);
    [~,hbrnew,suc]=C1branch_from_C2point(funcs,hh,hopf_wbifs.parameter.free,...
        'codim1','hopf','codim2','hoho',hopfopts{:});
    assert(all(suc(:)>0))
    indh=indh+1;
    [hbr(indh),succ,fail,rjct]=br_contn(funcs,hbrnew(1),nhsteps,'plotaxis',ax2);
    assert(succ>0&&rjct==0&&fail<succ/4)
    hbr(indh)=br_rvers(hbr(indh));
    [hbr(indh),succ,fail,rjct]=br_contn(funcs,hbr(indh),nhsteps,'plotaxis',ax2);
    assert(succ>0&&rjct==0&&fail<succ/4)
end
%% Branch off to Torus bifs at Hopf-Hopf interactions
ntrsteps=5;
trbr=repmat(triv_eqs,1,0);
for i=1:length(hohos)
    hh=hohos(i);
    [trfuncs,trbrnew,suc]=C1branch_from_C2point(funcs,hh,hopf_wbifs.parameter.free,...
        'codim1','TorusBifurcation','codim2','hoho','print_residual_info',1,hopfopts{:},'matrix','sparse');
    assert(all(suc(:)>0))
    [trbr(end+1),succ,fail,rjct]=br_contn(trfuncs,trbrnew(1),ntrsteps,'plotaxis',ax2);
    assert(succ>0)
    [trbr(end+1),succ,fail,rjct]=br_contn(trfuncs,trbrnew(2),ntrsteps,'plotaxis',ax2);
    assert(succ>0)
end
