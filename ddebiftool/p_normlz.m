function p=p_normlz(p)
%% Normalize selected fields of point (during |br_contn|)
% function normalized_p=p_normlz(p)
% INPUT:
%	p point to be normalized
% OUTPUT:
%	normalized_p normalized point
%
% (c) DDE-BIFTOOL v. 2.00, 21/10/2001
%
% $Id: p_normlz.m 367 2019-07-14 21:56:46Z jansieber $
%%
p=dde_apply({'dde_',p.kind,'_normlz'},p);
end
