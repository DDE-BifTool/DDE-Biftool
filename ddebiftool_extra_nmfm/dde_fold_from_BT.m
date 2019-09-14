function fold=dde_fold_from_BT(point,varargin)
%% convert Bogdanov-Takens point to fold
%
% $Id: dde_fold_from_BT.m 308 2018-10-28 15:08:12Z jansieber $
%%
fold=dde_fold_create('point',point,'v',point.nvec.q0);
if isfield(point,'stability')
    fold.stability=point.stability;
end
end