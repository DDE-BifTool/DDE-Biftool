%%  DDE-BIFTOOL demo 1 - Neuron: run all scripts of neuron example for testing
% This script runs all demos files in a single run for testing purposes.
% See <html/demo1.html> for the published version of the single scripts
%
% $Id: publishdemo.m 367 2019-07-14 21:56:46Z jansieber $
%
%% Description and load path
publish('demo1');
%% Definition of user functions
opts={'maxOutputLines',20};
publish('demo1_funcs',opts{:});
%% Symbolic definition of right-hand side and derivatives
publish('gen_sym_demo1');
%% Equilibria
publish('demo1_stst',opts{:});
%% Hopf bifurcations
publish('demo1_hopf',opts{:});
%% Normal forms along Hopf bifurcations (broken for r109 of br_bifdet)
% publish('demo1_normalforms','maxOutputLines',Inf);
%% Periodic orbits
publish('demo1_psol',opts{:});
%% Homoclinic connection
publish('demo1_hcli',opts{:});
%% Saddle-node of periodic orbits
publish('demo1_POfold',opts{:});
%% repeat procedure using convenience functions
publish('demo1_simple',opts{:});
%% Background information
publish('demo1_definingsystems');