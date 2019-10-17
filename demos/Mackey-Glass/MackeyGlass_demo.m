%% DDE-Biftool demo Mackey-Glass Equation
%
% The Mackey-Glass equation is given by
% 
% $$x'(t)=\beta \frac{x(t-\tau)}{1+x(t-\tau)^n}-\gamma x(t)$$
% 
% Parameters are (in this order) |beta|, |n|, |tau| (|gamma| is not part of
% parameter vector).
%
% <html>
% $Id: MackeyGlass_demo.m 178 2017-03-14 00:00:01Z jansieber $
% </html>
%
%% load DDE-Biftool into path
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');
format compact
%% Initial parameters and state
parnames={'beta','n','tau','gamma'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
gamma=1.0;
beta=2;
n=10;
tau=0;
x0=(beta-1)^(1/n);
par0=cellfun(@(x)evalin('caller',x),parnames); %[beta,n,tau,gamma];
bounds={'max_bound',[ind.tau,2;ind.beta,5],'max_step',[0,0.3]};
%% Set user-defined functions
%  Using |beta|,|n|,|tau| and |gamma| as parameters. Different ways to
%  define the right-hand side: |fvec| uses finite differences for
%  derivatives in a vectorized manner. |fsingle| does the same but does not
%  employ vectorization (should be slower), |fsymbolic| uses the right-hand
%  side and derivatives created via symbolic toolbox (see
%  <gen_sym_MackeyGlass.html>).
rhs=@(x,xtau,beta,n)beta*xtau./(1+xtau.^n)-gamma*x;
fvec=set_funcs(...
    'sys_rhs',@(xx,p)rhs(xx(1,1,:),xx(1,2,:),p(ind.beta),p(ind.n)),...
    'sys_tau',@()ind.tau,'x_vectorized',true);
fsingle=set_funcs(...
    'sys_rhs',@(xx,p)rhs(xx(1,1,:),xx(1,2,:),p(ind.beta),p(ind.n)),...
    'sys_tau',@()ind.tau);
fsymbolic=set_symfuncs(@sym_MackeyGlass_mf,'sys_tau',@()ind.tau);
%% Choose which problem definition should be tested
funcs=fsymbolic;
%% Initialization of branch of non-trivial equilibria
contpar=ind.tau;
nontriv_eqs=SetupStst(funcs,'x',x0,'parameter',par0,'step',0.1,...
    'contpar',contpar,'max_step',[contpar,0.3],bounds{:});
%% Compute, find stability and bifurcations of non-trivial equilibria 
disp('Trivial equilibria');
figure(1);clf;ax1=gca;
nontriv_eqs=br_contn(funcs,nontriv_eqs,3,'plotaxis',ax1);
[nontriv_eqs,~,ind_hopf,bif1types]=LocateSpecialPoints(funcs,nontriv_eqs);
nunst_eqs=GetStability(nontriv_eqs);
fprintf('Hopf bifurcation near point %d\n',ind_hopf);
%% Continue Hopf bifurcation in two parameters
[hbranch,suc]=SetupHopf(funcs,nontriv_eqs,ind_hopf,...
    'contpar',[ind.beta,ind.tau],'dir',ind.beta,'step',1e-1,bounds{:});
figure(2);clf;ax2=gca;
hbranch=br_contn(funcs,hbranch,30,'plotaxis',ax2);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,30,'plotaxis',ax2);
%% Compute L1 coefficient 
% to find if Hopf bifurcation is supercritical (L1<0) or subcritical (L1>0)
[hbranch,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hbranch)
fprintf('maximal L1 coefficient along Hopf branch: %g\n',max(hopftests.genh(1,:)));
fprintf('max of error estimate for L1 coefficient: %g\n',norm(diff(hopftests.genh),'inf'));
%% Branch off at  Hopf bifurcation
disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');
[per_orb,suc]=SetupPsol(funcs,nontriv_eqs,ind_hopf,...
    'print_residual_info',1,'intervals',20,'degree',4,...
    'max_bound',[contpar,20],'max_step',[contpar,0.5],'matrix','full');
if ~suc
    error('MackeyGlassDemo:fail',...
        'MackeyGlassDemo: initialization of periodic orbit failed');
end
figure(1);
hold on
per_orb=br_contn(funcs,per_orb,60);
[nunst_per,dom_per,triv_per,per_orb.point]=...
    GetStability(per_orb,'exclude_trivial',true,'funcs',funcs);
%% Find period doubling bifurcations in two parameters
ind_pd=find(diff(nunst_per)==1);
[pdfuncs,pdbranch1,suc]=SetupPeriodDoubling(funcs,per_orb,ind_pd(1),...
    'contpar',[ind.beta,ind.tau],'dir',ind.beta,'step',1e-1,bounds{:});
if ~suc
    error('MackeyGlassDemo:fail',...
        'MackeyGlassDemo: initialization of period doubling failed');
end
figure(2);
pdbranch1=br_contn(pdfuncs,pdbranch1,30);
pdbranch1=br_rvers(pdbranch1);
pdbranch1=br_contn(pdfuncs,pdbranch1,30);
%% Check Floquet multipliers 
% (note that Floquet multipliers are often unreliable)
[nunst_pd,floqpd1,triv_defect1,pdbranch1.point]=GetStability(pdbranch1,...
    'exclude_trivial',true,'funcs',pdfuncs); 
fprintf('max defect of Floquet multiplier at -1 or 1: %g\n',max(abs(triv_defect1)));
%% Branch off at period doubling 
% (Solutions at far end get inaccurate.)
[per2,suc]=DoublePsol(funcs,per_orb,ind_pd(1));
if ~suc
    error('MackeyGlassDemo:fail',...
        'MackeyGlassDemo: branching off at period doubling failed');
end
figure(1);
per2=br_contn(funcs,per2,60);
[nunst_per2,dom,triv_defect]=GetStability(per2,'funcs',funcs,'exclude_trivial',true); 
fprintf('max defect of Floquet multiplier at 1: %g\n',max(triv_defect));
%% Continue period doublings in two parameters for secondary PD
ind_pd2=find(diff(nunst_per2)==1);
[pd2funcs,pdbranch2,suc]=SetupPeriodDoubling(funcs,per2,ind_pd2(1),...
    'contpar',[ind.beta,ind.tau],'dir',ind.beta,'step',1e-1,bounds{:});
if ~suc
    error('MackeyGlassDemo:fail',...
        'MackeyGlassDemo: initialization of 2nd period doubling failed');
end
figure(2);
pdbranch2=br_contn(pdfuncs,pdbranch2,30);
pdbranch2=br_rvers(pdbranch2);
pdbranch2=br_contn(pdfuncs,pdbranch2,30);
%% Check Floquet multipliers along period doubling bifurcation
% (Note that Floquet multipliers are often unreliable.)
[nunst_pd2,floqpd2,triv_defect2,pdbranch2.point]=GetStability(pdbranch2,...
    'exclude_trivial',true,'funcs',pdfuncs);
fprintf('max defect of Floquet multiplier at -1 or 1: %g\n',max(abs(triv_defect2)));

%% Two-parameter bifurcation diagram
% Assigning a name and color to a curve. Others are chosen automatically
figure(3)
clf;ax3=gca;hold(ax3,'on')
lg=Plot2dBranch(hbranch);
lg=Plot2dBranch(pdbranch1,'lgname','PD1','funcs',pdfuncs,'oldlegend',lg);
lg=Plot2dBranch(pdbranch2,'lgname','PD2','funcs',pdfuncs,'oldlegend',lg,'color',[0,0.5,0]);
xlabel('\beta');
ylabel('\tau');
title(sprintf(['Bifurcation diagram of Mackey-Glass eqn,\nother parameters, ',...
    ' n=%g, gamma=%g'],par0(ind.n),par0(ind.gamma)));
grid on
%% Time profiles of period doubling orbits and errors of Floquet multipliers
bifsols={pdbranch1.point,pdbranch2.point};
floqpd={triv_defect1,triv_defect2};
get_par=@(i,k)arrayfun(@(x)x.parameter(i),bifsols{k});
figure(4)
clf;
for k=1:2
    subplot(2,2,k);
    hold on
    for i=1:length(bifsols{k})
        plot(bifsols{k}(i).mesh*bifsols{k}(i).period,bifsols{k}(i).profile(1,:),'-');
    end
    hold off
    box on
    grid on
    title(sprintf('PD%d: time profiles of period doubling',k));
    xlabel('t');
    ylabel('x');
    subplot(2,2,2+k);
    semilogy(1:length(bifsols{k}),floqpd{k},'.-');
    grid on
    title(sprintf('PD%d: dist crit Floq mult to 1 or -1',k));
    ylabel('error');
    xlabel('point along branch');
end
%% save
save('MGresults.mat');
