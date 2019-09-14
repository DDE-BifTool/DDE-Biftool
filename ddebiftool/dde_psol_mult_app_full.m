%% construct eigenvalue problem for Floquet multipliers using full matrices
%
% $Id: dde_psol_mult_app_full.m 351 2019-06-19 14:38:12Z jansieber $
%%
function [s,ef]=dde_psol_mult_app_full(Marg,varargin)
default={'geteigenfuncs',false,'method',[]};
options=dde_set_options(default,varargin,'pass_on','method');
ef=[];
if ~options.geteigenfuncs
    s=eig(Marg{:});
else
    [ef,s]=eig(Marg{:});
    s=diag(s);
end
end
