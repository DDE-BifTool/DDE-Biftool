%% continue Hopf bifurcations
%
%
clear
load('humphries1dbifeqs');
indhopf=find(eqnunst>0,1,'first');
[hbranch1,suc]=SetupHopf(funcs,eqbr,indhopf,'contpar',[1,2],'dir',1,'step',0.1);
if ~suc
    error('Hopf initialization failed');
end
clf
hbranch1=br_contn(funcs,hbranch1,20);
hbranch1=br_rvers(hbranch1);
hbranch1=br_contn(funcs,hbranch1,60);
%% continue equilibrium in kappa2 to find other dominant Hopf bifurcation
p0(1)=1;
[eqbr2,suc]=gen_stst(funcs,'contpar',2,'x',0,'parameter',p0,...
    'max_bound',[2,12],'max_step',[2,0.1]);
clf
eqbr2=br_contn(funcs,eqbr2,50);
%% stability
eqbr2.method.stability.minimal_real_part=-5;
[eq2nunst,dom,triv_defect,eqbr2.point]=...
GetStability(eqbr2,'funcs',funcs,'points',1:length(eqbr2.point));
%%
indhopf2=find(eq2nunst>0,1,'first');
[hbranch2,suc]=SetupHopf(funcs,eqbr2,indhopf2,'contpar',[1,2],'dir',1,'step',0.1);
if ~suc
    error('Hopf initialization failed');
end
clf
hbranch2=br_contn(funcs,hbranch2,40);
hbranch2=br_rvers(hbranch2);
hbranch2=br_contn(funcs,hbranch2,10);
%%
save('hopfs.mat')
