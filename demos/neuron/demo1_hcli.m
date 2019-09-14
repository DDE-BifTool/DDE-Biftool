%% Continuation of homoclinic connections
%
% <html>
% $Id: demo1_hcli.m 367 2019-07-14 21:56:46Z jansieber $
% </html>
%
% (requires running <demo1_psol.html> first to create branch5)
%%
%#ok<*ASGLU,*NOPTS,*NASGU>
%
%% Cutting a long-period periodic orbit
% The 7th last point along the branch of periodic orbits |branch5| was
% close to a homoclinic connection, consisting of two loops. Its period is
% also not too long such that the coarse mesh gives still an accurate
% solution. Using the (added) routines to compute homoclinic solutions, we
% correct each of the two loops to a homoclinic orbit, thereby obtaining
% also some stability information of the steady state point. We take the
% first half of the profile and rescale it to $[0,1]$. Then we convert it
% into a homoclinic/heteroclinic structure (point), and correct it with a
% Newton iteration (|p_correc|).
hcli1=psol;
hcli1.mesh=hcli1.mesh(1:65);
hcli1.profile=hcli1.profile(:,1:65);
hcli1.period=hcli1.period*hcli1.mesh(end);
hcli1.mesh=hcli1.mesh/hcli1.mesh(end);

hcli1=p_tohcli(funcs,hcli1) % convert it to a point of homoclinic structure

mh=df_mthod(funcs,'hcli');
[hcli1,success]=p_correc(funcs,hcli1,ind_a21,[],mh.point) % and correct it

%% The second half of the profile
% We apply the same procedure on the second half of the profile.
hcli2=psol;
hcli2.mesh=hcli2.mesh(81:end-16);
hcli2.profile=hcli2.profile(:,81:end-16);
hcli2.mesh=hcli2.mesh-hcli2.mesh(1);
hcli2.period=hcli2.period*hcli2.mesh(end);
hcli2.mesh=hcli2.mesh/hcli2.mesh(end);

hcli2=p_tohcli(funcs,hcli2);
[hcli2,success]=p_correc(funcs,hcli2,ind_a21,[],mh.point);
figure(14);clf;subplot(2,1,1);
p_pplot(hcli1);
xlabel('t/period');ylabel('x1, x2');
subplot(2,1,2);
p_pplot(hcli2);
xlabel('t/period');ylabel('x1, x2');
%% Figure: time profiles of both connection loop separately
% Homoclinic profiles of the two loops depicted in figure
% <demo1_psol.html#longperiod>, now computed using the defining system for
% homoclinic/heteroclinic connections.
%% Mesh Refinement
% We recompute the first homoclinic orbit, using 70 intervals, and correct
% this point.
hcli1=p_remesh(hcli1,4,70);
[hcli1,success]=p_correc(funcs,hcli1,ind_a21,[],mh.point)
%% Continuation of homoclinic in two parameters
% If we free a second parameter, we can continue this homoclinic orbit with
% respect to two free parameters.  As a second free parameter, we choose
% $\tau_s$. We first create a default branch of homoclinic orbits, add
% |hcli1| as a first point, perturb it, and add the corrected
% perturbation as a second point.
figure(15);
branch6=df_brnch(funcs,[ind_a21 ind_taus],'hcli');
branch6.point=hcli1;
hcli1.parameter(ind_taus)=1.49;
[hcli1,success]=p_correc(funcs,hcli1,ind_a21,[],mh.point);
branch6.point(2)=hcli1;
[branch6,s,r,f]=br_contn(funcs,branch6,19);
xlabel('a21');ylabel('\tau_s');
%% Figure: Two parameter bifurcation diagram with homoclinic connection
save('demo1_hcli_results.mat');
