%% Cusp demo - generate right-hand side and derivatives from symbolic expressions
%
% <html>
% $Id$
% </html>
%
% Demo contributed by Maikel Bosschaert, extended and modified by Jan
% Sieber
% 
% from Giannakopoulos, F. and Zapp, A. (2001). Bifurcations in a planar
% system of differential delay equations modeling neural activity. Physica
% D: Nonlinear Phenomena, 159(3):215-232.
%
%% Differential equations
%
% $$\mu\dot{u}_{1}(t)=-u_{1}(t)+q_{11}\alpha(u_{1}(t-T))-q_{12}u_{2}(t-T)+e_{1}$$
%
% $$\mu\dot{u}_{2}(t)=-u_{2}(t)+q_{21}\alpha(u_{1}(t-T))-q_{22}u_{2}(t-T)+e_{2}$$
%% include path to necessary functions
clear
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_symbolic/');
%% Define symbolic variables in their dimensions
ntau=1;
xx=sym('u',[2,ntau+1]);
par=sym('p',[1,6]);

alpha = 1./(1+exp(-4*xx(1,2)))-1/2;

f=[ -xx(1,1)+par(1)*alpha-par(2)*xx(2,2)+par(4);...
    -xx(2,1)+par(3)*alpha+par(5)];
%% generate right-hand side and directional derivatives
% all functions will be stored in file with name sym_cusp_rhs.m
dde_sym2funcs(f,xx,par,'filename','sym_cusp');