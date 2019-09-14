function [stst,stpcnd]=p_tostst(point,varargin)
%% convert point to stst (from fold, hopf, hcli or near branching)
% function [stst_point,stpcnd]=p_tostst(funcs,point,pert)
% INPUT:
%	point 
%       % for branch switching from stst
%       optional:
%       'radius': size of the perturbation onto the emanating branch
%       'funcs' problem functions

% OUTPUT:
%	stst_point starting guess for stst point derived from point
% COMMENT:
%       if point kind is stst and the point is near a 
%         branching bifurcation, stst_point is a guess on the
%         branch and stpcnd should be used as an extra condition
%         in p_correc 
%       if point kind is fold or hopf, stst_point contains the
%         same point, but now as steady state only, no stpcnd
%         is returned
%       if point kind is hcli stst_point(1) is initial
%         and stst_point(2) is final steady state of the
%         given connecting orbit, no stpcnd is returned

% (c) DDE-BIFTOOL v. 2.03, 01/08/2002
%
% $Id: p_tostst.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
%% for backward compatibility re-organize arguments
args=varargin;
if ~isfield(point,'kind')
    funcs=point;
    point=args{1};
    args=[args(2:end),{'funcs',funcs}];
end
switch point.kind
    case 'stst'
        [stst,stpcnd]=dde_stst_from_stst(point,args{:});
    case 'hcli'
        stst(2)=dde_stst_create('point',point,'x',point.x2);
        stst(1)=stst(2);
        stst(1).x=point.x1;
    otherwise
        try
            stst=dde_apply({'dde_stst_from_',point.kind,''},point,args{:});
        catch
            stst=dde_stst_create('point',point);
        end
end
end
