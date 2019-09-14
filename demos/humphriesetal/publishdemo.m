%%  DDE-BIFTOOL demo for sd-DDE: run all scripts of example from Humphries etal for testing
% This script runs all demos files in a single run for testing purposes.
% See <html/humphriesetal_demo.html> for the published version of the
% single scripts
%
% <html>
% $Id: publishdemo.m 176 2017-03-13 00:25:33Z jansieber $
% </html>
%
%% Definition of user functions and load path
publish('humphriesetal_demo');
%%
%% Equilibria and their bifurcations with normal forms
publish('humphriesetal_equilibria');
%% Periodic orbits
opts={'maxOutputLines',20};
publish('humphriesetal_periodic1dbif',opts{:});
%% Bifurcations of periodic orbits
opts={'maxOutputLines',20};
publish('humphriesetal_periodic2dbif',opts{:});
