%% One-parameter periodic orbit continuation of example in Humphries et al (DCDS-A 2012)
%
% <html>
% $Id$
% </html>
%
%%  Periodic orbits branching off from 1st Hopf bifurcation
% We continue the periodic orbits in |parameter(1)| (|kappa1|), and calculate
% their stability.
eqnunst=GetStability(eqbr);
indhopf=find(eqnunst>0,1,'first');
[per,suc]=SetupPsol(funcs,eqbr,indhopf,'contpar',ind.kappa1,'degree',5,'intervals',30,...
    'print_residual_info',1,'radius',1e-2);
if ~suc
    error('initialization of periodic orbits failed');
end
per.parameter.max_step=[ind.kappa1,0.3];
figure(1);ax1=gca;hold(ax1,'on');
per=br_contn(funcs,per,200,'plotaxis',ax1);
disp('calculate stability')
[pernunst,dom,triv_defect,per.point]=...
    GetStability(per,'exclude_trivial',true,'funcs',funcs); %#ok<*ASGLU>
fprintf('maximum error of trivial Floquet multiplier: %g\n',max(abs(triv_defect)));
%% One-parameter bifurcation diagram for family of periodic orbits
% Continuation parameter is $\kappa_1$. Add to diagram for  equilibria.
figure(3);ax3=gca;hold(ax3,'on');
Plot2dBranch(per);
grid on
xlabel('\kappa_1');
ylabel('amplitude');
lg=legend(ax3);
set(lg,'location','northwest');
axis(ax3,'auto')
grid on
%% Next step: continuation of bifurcations of periodic orbits in a two parameters
% see <humphriesetal_periodic2dbif.html>.
save('humphriesetal_1dbif.mat');
