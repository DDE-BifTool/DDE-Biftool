function fold=dde_fold_from_zeho(zeho,varargin)
%% wrapper to extract fold point from zeho
%
% $Id: dde_fold_from_zeho.m 308 2018-10-28 15:08:12Z jansieber $
%%
if isempty(zeho.nvec)
    fold=dde_fold_from_stst(zeho,varargin{:});
else
    fold=dde_fold_create('point',zeho,'v',zeho.nvec.q0,varargin{:});
end
end

