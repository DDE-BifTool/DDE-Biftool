%% first Lyapunov coefficient L1 in Hopf point along Hopf branch
% (wrapper for NormalformCoefficients)
%% Input
%
% * funcs: problem functions
% * branch: hopf branch or array of points of kind 'hopf'
% * optional name-value pairs: 
% * 'use_nullvectors' (boolean, default false) use nullvectors from
% previous point along branch to compute new nullvectors by bordering
% (otherwise, svd is used)
%
%% Output
%
% * L1: array of L1 coefficients along branch
% * L1low: equals L1 if funcs.sys_mfderi is provided, otherwise (for
% finite-difference approximation) a low-er-order estimate of L1
% coefficients. Use L1-L1low to estimate the accuracy of L1. If this is
% large, finite-differences are problematic. If this is small it may(!) be
% ok.
%%
function [L1,L1low,branch]=HopfLyapunovCoefficients(funcs,branch,varargin)
%%
%
% $Id: HopfLyapunovCoefficients.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
[L1,L1low,branch]=NormalformCoefficients(funcs,branch,varargin{:});
end