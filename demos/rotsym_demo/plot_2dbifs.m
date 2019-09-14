%% plot 2d bifurcation diagram
%
% $Id: plot_2dbifs.m 354 2019-06-30 23:15:16Z jansieber $
%
clear
load('LKbifs');
figure(4);clf
ppd=getpar(ip.phi,pdbr);
epd=getpar(ip.eta,pdbr);
ppf=getpar(ip.phi,mwfoldbr);
epf=getpar(ip.eta,mwfoldbr);
plot([-4,8],par(ip.eta)*[1,1],'k-','linewidth',2);
hold on
plot(ph,eh,'ro-',...
    ppd,epd,'ks-','linewidth',2);
plot(ppf,epf,'d-','color',[0,0.5,0],'linewidth',2);
plot( pf,ef,'b.-', pf2,ef2,'b.-', pf2+2*pi,ef2,'b.-','linewidth',2);
axis([-4,8,0,0.01]);
grid on
set(gca,tdeco{:});
legend({'1d bif','RWHopf','MW PD','MW,fold','RWfold'},'location','eastoutside');
xlabel('$\phi$',ldeco{:});
ylabel('$\eta$',ldeco{:});
title('bifs of rotating and modulated waves in Lang-Kobayashi system',tdeco{:});
