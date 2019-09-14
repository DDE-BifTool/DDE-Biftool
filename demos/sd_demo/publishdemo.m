%% DDE-BIFTOOL  state-dependent delays sd-demo
% This script publishes all demos files in a single run for testing purposes.
% See <html/sd_demo.html> for the published version of the single scripts
%
% $Id: publishdemo.m 20 2014-04-11 19:27:33Z jan.sieber $
%
%% Description and load path, definition of user functions, bifurcation analysis
publish('sd_demo.m');
publish('gen_sym_sd_demo.m','evalCode',false);
