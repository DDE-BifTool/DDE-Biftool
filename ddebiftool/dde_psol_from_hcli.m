function [psol,tangent]=dde_psol_from_hcli(point,varargin)
%% convert connecting orbit to psol, assuming long period
% tangent assumes change only in period
%
% $Id: dde_psol_from_hcli.m 308 2018-10-28 15:08:12Z jansieber $
%%
psol=dde_psol_create('point',point,varargin{:});
tangent=dde_psol_create('point',psol,'profile',0*psol.profile,'period',1);
end
