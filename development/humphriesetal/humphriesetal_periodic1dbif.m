%% One-parameter periodic orbit continuation of example in Humphries et al (DCDS-A 2012)
%
% <html>
% $Id$
% </html>
%
%%  Periodic orbits branching off from 1st Hopf bifurcation
% We continue the periodic orbits in |parameter(1)| (|kappa1|), and calculate
% their stability.
indhopf=find(eqnunst>0,1,'first');
[per,suc]=SetupPsol(funcs,eqbr,indhopf,'contpar',indkappa1,'degree',5,'intervals',30,...
    'print_residual_info',1,'radius',1e-2);
if ~suc
    error('initialization of periodic orbits failed');
end
per.parameter.max_step=[indkappa1,0.3];
per=br_contn(funcs,per,200);
disp('calculate stability')
[pernunst,dom,triv_defect,per.point]=...
    GetStability(per,'exclude_trivial',true,'funcs',funcs); %#ok<*ASGLU>
fprintf('maximum error of trivial Floquet multiplier: %g\n',max(abs(triv_defect)));
%% One-parameter bifurcation diagram for family of periodic orbits
% Continuation parameter is $\kappa_1$.
ppars=arrayfun(@(x)x.parameter(1),per.point);
pmeshes=cell2mat(arrayfun(@(x)x.mesh(:),per.point,'uniformoutput',false));
pprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',per.point,'uniformoutput',false));
clf
clrs=colormap('lines');
pernunst_cases=unique(pernunst);
amp=max(pprofs)-min(pprofs);
hold on
lstr={};
for i=1:length(pernunst_cases);
    sel=pernunst==pernunst_cases(i);
    plot(ppars(sel),amp(sel),'o','color',clrs(i,:));
    lstr={lstr{:},sprintf('#unst=%d',pernunst_cases(i))}; %#ok<CCAT>
end
hold off
grid on
xlabel('\kappa_1');
ylabel('amplitude');
legend(lstr,'location','southeast');
grid on
%%
save('humphries1dbif.mat');
