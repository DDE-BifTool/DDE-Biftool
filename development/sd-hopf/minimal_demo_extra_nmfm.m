%% Minimal demo - Normal forms of Hopf bifurcations
% This part creates the computations of normal form coefficients along Hopf
% bifurcations, requiring theextension |ddebiftool_extra_nmfm|. This demo
% requires  <minimal_demo_stst_psol.html> to have
% run beforehand.
%
% <html>
% $Id: minimal_demo_extra_nmfm.m 109 2015-08-31 23:45:11Z jansieber $
% </html>
%
%%
%#ok<*SAGROW>
load('minimal_demo_stst_psol_results.mat');
hopfcurves={hopf,hopf1};
ind_hoho=0;
for i=1:length(hopfcurves)
    %% Compute Lyapunov coefficient L1 along Hopf curves
    [L1{i},L1low]=HopfLyapunovCoefficients(funcs,hopfcurves{i});
    fprintf('maximal error of L1 along hopf branch=%g\n',norm(L1{i}-L1low,'inf'));
    %% Hopf bifurcation changes criticality
    % Detect generalized Hopf bifurcation
    [genh{i},genhlow,hopfref{i},ind_genh(i)]=GeneralizedHopfNormalform(funcs,hopfcurves{i},...
        find(diff(sign(L1{i}))~=0,1)+(0:1));
    fprintf(['Generalized Hopf point at (b,tau)=(%g,%g)\n',...
        'with L1=%g, L2=%g,\n',...
        'L1 error est=%g, L2 error est=%g.\n'],...
        genh{i}.parameter(indb),genh{i}.parameter(indtau),...
        genh{i}.nmfm.L1,genh{i}.nmfm.L2,...
        abs(genh{i}.nmfm.L1-genhlow.nmfm.L1),abs(genh{i}.nmfm.L2-genhlow.nmfm.L2));
    %% Compute stability along Hopf curve
    % This shows that there are several HopfHopf interactions: detect them
    % and compute their normal form.
    [nunsth{i},dum,dum,hopfref{i}.point]=GetStability(hopfref{i},'funcs',funcs,...
        'exclude_trivial',true,'locate_trivial',@(p)[-1i*p.omega,1i*p.omega]); %#ok<ASGLU>
    ind_hh=find(abs(diff(nunsth{i}))==2);
    for k=1:length(ind_hh)
        ind_hoho=ind_hoho+1;
        [hoho{ind_hoho},hoho_low]=HopfHopfNormalform(funcs,hopfref{i},ind_hh(k)+(0:1)); 
        fprintf(['Normal form coefficients of Hopf-Hopf point\n',...
            'at (b,tau)=(%g,%g) with omega1=%g, omega2=%g:\n'],...
            hoho{ind_hoho}.parameter(indb),hoho{ind_hoho}.parameter(indtau),...
            hoho{ind_hoho}.omega1,hoho{ind_hoho}.omega2);
        disp(hoho{ind_hoho}.nmfm);
        fprintf('Error of normal form coefficients: %g\n',...
        norm(structfun(@(x)x,hoho{ind_hoho}.nmfm)-structfun(@(x)x,hoho_low.nmfm),'inf'));
    end
end
%% Save and continue
% For continuation of folds and torus bifurcations of periodic orbits, see
% <minimal_demo_extra_psol.html>. Final results in
% <minimal_demo_plot_2dbif.html>.
save('minimal_demo_extra_nmfm_results.mat')