%% plot & decorate 2d bifurcation diagram
clear
load('hopfs','hbranch1','hbranch2');
load('humphries_2dbif','trbr','pbr');
%%
getpar=@(br)cell2mat(arrayfun(@(x)x.parameter(1:2)',br.point,'uniformoutput',false));
tr=getpar(trbr);
hp1=getpar(hbranch1);
hp2=getpar(hbranch2);
pf=getpar(pbr);
clf
hold on
ldeco={'linewidth',2};
fdeco={'fontsize',16,'fontweight','bold'};
plot(hp1(1,:),hp1(2,:),'r',ldeco{:});
plot(hp2(1,:),hp2(2,:),'--','color',[1,1,1]*0.5,ldeco{:});
plot(pf(1,:),pf(2,:),'b-',ldeco{:});
plot(tr(1,:),tr(2,:),'-','color',[0,0.5,0],ldeco{:});
axis([0,12,0,5]);
legend({'Hopf1','Hopf2','PO Folds','Torus bif'},'location','best');
set(gca,'box','on',ldeco{:},fdeco{:});
print('-depsc2','-r300','humphriesetal.eps');
