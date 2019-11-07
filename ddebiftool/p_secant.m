function secant=p_secant(secant,norm_point)
%% compute normalized secant
% function n_secant=p_secant(secant,norm_point)
% INPUT: 
%	secant unnormalized secant
%	norm_point norm of a computed reference point
% OUTPUT:
%	n_secant normalized secant

% (c) DDE-BIFTOOL v. 2.00, 21/10/2001 
%
% $Id: p_secant.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
secant=p_axpy(norm_point/p_norm(secant),secant,[]);
end