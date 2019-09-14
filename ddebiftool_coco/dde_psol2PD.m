function data=dde_psol2PD(funcs,info,varargin)
%% constructor for period doubling bifurcation from psol
% wrapper for dde_psol2torus
% $Id: dde_psol2PD.m 346 2019-05-13 05:41:50Z jansieber $
data=dde_psol2torus(funcs,info,varargin{:},'biftype','PD','closest',-1);
end
