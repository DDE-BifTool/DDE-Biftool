function ip = p_inprod(p1, p2)
%% Inner product between p1 and p2
% function ip = p_inprod(p1, p2)
% PURPOSE: compute the inner product between p1 and p2
% INPUT:
%	p1, p2: points 
% OUTPUT:
%   ip: inner product of the points
%
% BW: the structure of this function is based on p_norm.
% The expressions were chosen such that p_inprod(p1,p1) = p_norm(p1)^2.
%
% $Id: p_inprod.m 115 2015-09-02 03:42:31Z jansieber $
%
%%
kind = p1.kind;
if ~strcmp(kind, p2.kind)
    % Not the same kind!
    error('P_INPROD: input points must be of the same kind!');
end

switch kind
    case 'stst'
        ip = [p1.parameter, p1.x']*[p2.parameter'; p2.x];
    case 'fold'
        ip = [p1.parameter, p1.x', p1.v']*[p2.parameter'; p2.x; p2.v];
    case {'hopf','genh','zeho','hoho'} % BW: Addition of genh and p1.v -> p1.v', p2.v' -> p2.v
        ip = [p1.parameter, p1.x', p1.v', p1.omega]*[p2.parameter'; p2.x; p2.v; p2.omega'];
    case 'psol'
        v1 = [p1.parameter, p1.profile/sqrt(size(p1.profile)), p1.period];
        v2 = [p2.parameter'; p2.profile'/sqrt(size(p2.profile)); p2.period];
        ip = v1*v2;
    case 'hcli'
        v1 = [p1.parameter, p1.profile/sqrt(size(size(p1.profile,2))), ...
            p1.period, p1.x1', p1.x2', p1.lambda_v, p1.lambda_w, p1.v/sqrt(size(p1.v,2)), p1.alpha];
        v2 = [p2.parameter'; p2.profile'/sqrt(size(size(p2.profile,2))); ...
            p2.period'; p2.x1; p2.x2; p2.lambda_v'; p2.lambda_w'; p2.v'/sqrt(size(p2.v,2)); p2.alpha];
        ip = v1*v2;
    otherwise
        display(kind);
        error('P_INPROD: point type is not recognized.');
end

if length(ip) ~= 1
    display(ip);
    error('P_INPROD: inner product does not result in scalar.');
end
end
