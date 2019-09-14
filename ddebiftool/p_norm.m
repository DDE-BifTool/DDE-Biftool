function n=p_norm(point,varargin)
%% norm of point for distance measuring
% function n=p_norm(point,...)
% INPUT:
%	point 
% OUTPUT:
%	n norm of point (equals sqrt(p_dot(point,point,...)))

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% $Id: p_norm.m 362 2019-07-14 15:49:40Z jansieber $
% 
%%
n=sqrt(p_dot(point,point,varargin{:}));
end
