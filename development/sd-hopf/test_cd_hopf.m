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
clear
load('minimal_demo_stst_psol_results.mat','hopf','funcs');
%% Compute Lyapunov coefficient L1 at one hopf point
funcs.x_vectorized=false;
tic;
hpl1=nmfm_func_hopf(funcs,hopf.point(1));
t1=toc
hpl2=nmfm_hopf(funcs,hopf.point(1));
t2=toc