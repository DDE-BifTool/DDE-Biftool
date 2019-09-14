%% SetupRWHopf - Initialize continuation of Hopf bifurcations of relative equilibria
%%
function [hbranch,suc]=SetupRWHopf(funcs,branch,ind,varargin)
%% 
% simple wrapper to have a sensible name (see demo rotsym_demo how to use this function).
%
% $Id: SetupRWHopf.m 309 2018-10-28 19:02:42Z jansieber $
%
[hbranch,suc]=SetupHopf(funcs,branch,ind,varargin{:});
end
