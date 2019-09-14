function [psol,tangent]=dde_psol_from_stst(point,varargin)
%% Create small amplitude harmonic oscillation around stst point after conversion to Hopf
%
% Output is initial periodic orbit and approximate tangent
% $Id: dde_psol_from_stst.m 308 2018-10-28 15:08:12Z jansieber $
%%
hopf=dde_hopf_from_stst(point,varargin{:});
[psol,tangent]=dde_psol_from_hopf(hopf,varargin{:});
end
