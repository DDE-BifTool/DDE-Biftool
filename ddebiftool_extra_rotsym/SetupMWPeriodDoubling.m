%% SetupMWPeriodDoubling - Initialize continuation of period doubling bifurcations of relative periodic orbits
%%
function [pdfuncs,pdbranch,suc]=SetupMWPeriodDoubling(funcs,branch,ind,varargin)
%%
% simple wrapper to have a sensible name and set number of trivial Floquet
% multipliers to two (see demo rotsym_demo how to do this).
%
% $Id: SetupMWPeriodDoubling.m 360 2019-07-04 23:11:26Z jansieber $
%
[pdfuncs,pdbranch,suc]=SetupMWTorusBifurcation(funcs,branch,ind,'biftype','PD',varargin{:});
end
