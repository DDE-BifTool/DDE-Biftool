%% Bifurcations of periodic orbits in two parameters for example from Humpries etal (DCDS-A 2012)
%
% <html>
% $Id$
% </html>
% 
%% Fold of periodic orbits
% Continuation parameters are $\kappa_1$ and $\kappa_2$. We remove stepsize
% restrictions for parameters, increase the maximum number of Newton
% iterations. The stability along the fold may indicate codimensions.
% However, increasing |triv_defect| indicates increasing errors in Floquet
% multiplier computations.
pf_ind0=find(diff(pernunst)==1,1,'first')+1;
per.method.point.print_residual_info=1;
per.parameter.max_step=[];
per.method.point.newton_max_iterations=8;
[pfuncs,pbr,suc]=SetupPOfold(funcs,per,pf_ind0,'contpar',[ind.kappa1,ind.kappa2],...
    'dir',ind.kappa2,'step',-0.01);
if ~suc
    error('initialization of fold of periodic orbits failed');
end
%%
figure(2);ax2=gca;
pbr=br_contn(pfuncs,pbr,220,'plotaxis',ax2);
pbr=br_rvers(pbr);
pbr=br_contn(pfuncs,pbr,60,'plotaxis',ax2);
[pfstab,dom,triv_defect,pbr.point]=GetStability(pbr,'exclude_trivial',true,'funcs',pfuncs);
pfpars=cell2mat(arrayfun(@(x)x.parameter(1:2)',pbr.point,'uniformoutput',false));
pfmeshes=cell2mat(arrayfun(@(x)x.mesh(:),pbr.point,'uniformoutput',false));
pfprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',pbr.point,'uniformoutput',false));
pf_amp=max(pfprofs)-min(pfprofs);
save('humphries_pofold.mat');
%% Torus bifurcation continuation
% Continuation parameters are $\kappa_1$ and $\kappa_2$. 
tr_ind0=find(diff(pernunst)==2,1,'first');
[trfuncs,trbr,suc]=SetupTorusBifurcation(funcs,per,tr_ind0,...
    'contpar',[ind.kappa1,ind.kappa2],'dir',ind.kappa2,'step',-0.01);
if ~suc
    error('initialization of torus bifurcation failed');
end
%%
figure(2);
trbr=br_contn(trfuncs,trbr,80);
trbr=br_rvers(trbr);
trbr=br_contn(trfuncs,trbr,40);
%% Stability of periodic orbits along torus bifurcation
% The code below demonstrates how one can use |GetStability| and its optional
% argument |'locate_trivial'| to exclude the known critical Floquet
% multipliers from the stability consideration. Thus, the stability changes
% along the torus bifurcation help detect codimension-two bifurcations.
[trstab,dom,triv_defect,trbr.point]=GetStability(trbr,...
    'exclude_trivial',true,'funcs',trfuncs);
%% Extract parameters, meshes and profiles of orbits at torus bifurcation
trpars=cell2mat(arrayfun(@(x)x.parameter(1:2)',trbr.point,'uniformoutput',false));
trmeshes=cell2mat(arrayfun(@(x)x.mesh(:),trbr.point,'uniformoutput',false));
trprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',trbr.point,'uniformoutput',false));
tr_amp=max(trprofs)-min(trprofs);
trmu=cell2mat(arrayfun(@(x)x.stability.mu(1:10),trbr.point,'uniformoutput',false));
%% 
% plot profiles of orbits at torus bifurcation
figure(5);clf
plot(trmeshes,trprofs);
grid on
xlabel('t/T');
ylabel('x');
%% Bifurcation diagram
figure(4);ax4=gca;hold(ax4,'on');
args={'ax',ax4,'stability',0.75,'markersize',4};
Plot2dBranch(pbr,'funcs',pfuncs,args{:});
Plot2dBranch(trbr,'funcs',trfuncs,args{:});
set(ax4,'ylim',[0,5]);
lg=legend(ax4);
set(lg,'location','eastoutside');
%%
save('humphriesetal_2dbif.mat');
