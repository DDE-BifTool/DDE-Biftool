%% SetupPeriodDoubling - Initialize continuation of period doubling bifurcation
%%
function [pdfuncs,pdbranch,suc]=SetupPeriodDoubling(funcs,branch,ind,varargin)
%% 
% Simple wrapper around SetupTorusBifurcation to have a sensible name and
% pass on bifurcation type to pdfuncs.
% See <SetupTorusBifurcation.html> for description of input and output.
%
% <html>
% $Id: SetupPeriodDoubling.m 309 2018-10-28 19:02:42Z jansieber $
% </html>
[pdfuncs,pdbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,'biftype','PD','closest',-1,varargin{:});
end
