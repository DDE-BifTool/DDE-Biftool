function fold=p_tofold(point,varargin)
%% convert stst point to fold point
% function fold_point=p_tofold(point,...)
% INPUT:
%	point point near fold point
%   optional 'funcs': problem functions
% OUTPUT:
%	fold_point starting guess for fold point derived from point

% (c) DDE-BIFTOOL v. 1.00, 08/04/2000
%
% $Id: p_tofold.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
%% for backward compatibility re-organize arguments
args=varargin;
if ~isfield(point,'x')
    funcs=point;
    point=args{1};
    args=[args(2:end),{'funcs',funcs}];
end
fold=dde_fold_from_stst(point,args{:});
end
