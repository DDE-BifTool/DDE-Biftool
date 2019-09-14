%% This file generates the rhs file, the standard derivative file,
%  the vectorized standard derivative file,
%  the multilinear forms needed for calculating 
%  the critical and parameter-dependent normal form coefficients
%  This script has been tested on matlab 2016b 
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: test_gen_mf.m 168 2017-03-03 22:04:32Z jansieber $
%
clear
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities/');
%% Define system
n = 3; % number of variables
taus = 1; % number of nonzero delays
p = 2; % number of parameters
free_pars = [1 2];

xx = sym('xx', [n,taus+1]);
par = sym('par', [1,9]);
v = sym('v', [1,n]);
taus = [0 par(3)];
m = length(taus);

% system
f=[...
    xx(2,1,:)-par(3)*xx(1,1,:).^3+par(4)*xx(1,2,:).^2-par(5)*xx(3,1,:)+par(1);
 par(5)-par(6)*xx(1,1,:).^2-xx(2,1,:);
 par(8)*(par(2)*(xx(1,1,:)-par(7))-xx(3,1,:))];

%%
[fstr,df,xxdev,pdev]=dde_sym2funcs(f,xx,par);