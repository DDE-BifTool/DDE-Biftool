%% Test state-dependent delay equations with three levels of nesting
%
% <html>
% $Id: nested_demo.m 20 2014-04-11 19:27:33Z jan.sieber $
% </html>
%
% The equation is 
% 
% $$x'(t)=-x(t-p_1-x(t-p_1-x(t-p_1-x(t))))+p_2x(t)^5$$
% 
% The parameter |p(1)| controls the delay at the Hopf bifurcation, |p(2)|
% controls stability of the periodic orbits at sufficiently large amplitude
% without influencing criticality of the Hopf bifurcation.
%%
clear
close all
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_psol/',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities/');
rmpath('../../ddebiftool_extra_nmfm/');
addpath('new_nmfm/');
%% Equation definition
% Change |ntau| to change the level of nesting.
for ntau=1
    funcs=set_funcs(...
        'sys_rhs',@(x,p)p-x(1,ntau+1,:),...
        'sys_ntau',@()ntau,...
        'sys_tau',@(nr,x,p)x(1,nr,:),'x_vectorized',true,...
        'sys_dirderi',@nested_dirderi);
    %% Trivial equilibrium
    % x=0 is always the only equlibrium.
    eqbr=SetupStst(funcs,'contpar',1,'x',pi/2,'parameter',pi/2-1e-2,'step',0.1);
    [eqnunst,~,~,eqbr.point]=GetStability(eqbr,'funcs',funcs);
    %% Determine Hopf bifurcation and compute criticality
    % At the Hopf bifurcation a family of periodic orbits branches off. Its
    % stability close to the Hopf bifurcation depends on the level of nesting |ntau|.
    hopf=SetupHopf(funcs,eqbr,1);
    hpl1=nmfm_hopf(funcs,hopf.point,'output',1);
    hpl2=nmfm_hopf(funcs,hopf.point,'output',2);
    L1{ntau}=hpl1.nmfm.L1
    L1low{ntau}=hpl2.nmfm.L1  
    per{ntau}=SetupPsol(funcs,eqbr,1,'plot',0);
    %% Periodic orbits continued in |p(1)|
    % The family of periodic orbits folds back and forth for |ntau=3|.
    per{ntau}.parameter.max_step=[0,0.01];
    per{ntau}=br_contn(funcs,per{ntau},10);
    %% Stability of periodic orbits
    % The only source of instability is the fold such that periodic orbits are
    % either stable or order-1 unstable.
    [pernunst{ntau},dom,triv_defect,per{ntau}.point]=...
        GetStability(per{ntau},'exclude_trivial',true,'funcs',funcs);
    fprintf('maximum error of trivial Floquet multiplier: %g\n',max(abs(triv_defect)));
end
%% Bifurcation diagrams periodic orbits
figure(1);clf;
hold on
lgtxt=cellfun(@(x,y)sprintf('k=%d, L1:%+8.5f',x,y),num2cell(1:ntau),L1,'uniformoutput',false);
clr=lines;
deco={'linewidth',2};
for i=1:ntau
    pmax=[pi/2,arrayfun(@(x)max(x.profile),per{i}.point)];
    pmin=[pi/2,arrayfun(@(x)min(x.profile),per{i}.point)];
    par=[pi/2,per{i}.point.parameter];
    pl{i}=plot(par,[pmin;pmax]','o-','color',clr(i,:),deco{:});
    lg(i)=pl{i}(1);
end
xlim=get(gca,'xlim');
plot(xlim,xlim,'k-',deco{:});
legend(lg,lgtxt,'location','best','fontname','fixedwidth','fontweight','bold','fontsize',12);
grid on
xlabel('p');
ylabel('max x, min x');
axis tight
set(gca,'box','on',deco{:},'fontname','fixedwidth','fontweight','bold','fontsize',12)
%% Save data
%print('-dpng','-r300','nestedhopf.png');
%save(sprintf('sd_basic_per%d.mat',ntau));