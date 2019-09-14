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
[pfuncs,pbr,suc]=SetupPOfold(funcs,per,pf_ind0,'contpar',[indkappa1,indkappa2],...
    'dir',indkappa2,'step',-0.01);
if ~suc
    error('initialization of fold of periodic orbits failed');
end
%%
figure(2);
pbr=br_contn(pfuncs,pbr,220);
pbr=br_rvers(pbr);
pbr=br_contn(pfuncs,pbr,60);
pforbits=pfuncs.get_comp(pbr.point,'solution');
[pfstab,dom,triv_defect,pforbitstab]=GetStability(pforbits,'exclude_trivial',true,...
    'locate_trivial',@(p)[1,1],'funcs',funcs);
pfpars=cell2mat(arrayfun(@(x)x.parameter(1:2)',pbr.point,'uniformoutput',false));
pfmeshes=cell2mat(arrayfun(@(x)x.mesh(:),pbr.point,'uniformoutput',false));
pfprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',pbr.point,'uniformoutput',false));
pf_amp=max(pfprofs)-min(pfprofs);
save('humphries_pofold.mat');
%% Torus bifurcation continuation
% Continuation parameters are $\kappa_1$ and $\kappa_2$. 
tr_ind0=find(diff(pernunst)==2,1,'first');
[trfuncs,trbr,suc]=SetupTorusBifurcation(funcs,per,tr_ind0,...
    'contpar',[indkappa1,indkappa2],'dir',indkappa2,'step',-0.01);
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
trorbits=trfuncs.get_comp(trbr.point,'solution');
trrot=trfuncs.get_comp(trbr.point,'omega');
trorbits=arrayfun(@(x,y)setfield(x,'parameter',[x.parameter,y]),trorbits,trrot);
trivial_floqs=@(p)[1,exp(1i*p.parameter(end)*pi),exp(-1i*p.parameter(end)*pi)];
[trstab,dom,triv_defect,trorbitstab]=GetStability(trorbits,'exclude_trivial',true,...
    'locate_trivial',trivial_floqs,'funcs',funcs);
%% 
% Extract parameters, meshes and profiles of orbits at torus bifurcation
trpars=cell2mat(arrayfun(@(x)x.parameter(1:2)',trbr.point,'uniformoutput',false));
trmeshes=cell2mat(arrayfun(@(x)x.mesh(:),trbr.point,'uniformoutput',false));
trprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',trbr.point,'uniformoutput',false));
tr_amp=max(trprofs)-min(trprofs);
trmu=cell2mat(arrayfun(@(x)x.stability.mu(1:10),trorbitstab,'uniformoutput',false));
%% 
% plot profiles of orbits at torus bifurcation
figure(3);clf
plot(trmeshes,trprofs);
grid on
xlabel('t/T');
ylabel('x');
%%
% bifucation diagram
figure(2);clf
hold on
unstabsel=trstab>0;
lw={'linewidth',2};
plot(trpars(1,unstabsel),trpars(2,unstabsel),'ro','markersize',4,'markerfacecolor','r',lw{:});
plot(hpars(1,:),hpars(2,:),'k-',...
    trpars(1,:),trpars(2,:),'r-',...
    pfpars(1,:),pfpars(2,:),'b-',lw{:});
set(gca,'xlim',[0,12],'ylim',[0,4]);
grid on
xlabel('\kappa_1');
ylabel('\kappa_2');
legend({'torus bif (unstab)','Hopf','torus bif','fold'},...
    'location','southwest');
%%
save('humphries_2dbif.mat');
