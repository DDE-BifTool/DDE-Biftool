function [y,J]=p_dot(point1,point2,varargin)
%% Compute dot product of two periodic solutions (psol structs)
% 2nd output J is derivative wrt to point1
% by default only the integral between profiles is taken
% point2 is remeshed if meshes are different (may make J invalid)
% optional inputs
%
% derivatives (2d-integer array [n1,n2], default [0,0]): compute scalar
%   product between n1th and n2th derivative of point1 and point2
% ind_comp (array of integers between 1 and size(point1.profile,1), default
%   is 1:size(point1.profile,1)): only include ind_comp components into
%   scalar product
% free_par_ind (array of integers between 1 and length(point1.parameter),
%   default is []): add
%   point1.parameter(free_par_ind)'*point2.parameter(free_par_ind) to scalar
%   product
% period (logical, default is false): include point1.period*point2.period
% into scalar product
%
% $Id: p_dot.m 348 2019-06-19 13:09:01Z jansieber $
%

%% process options
[y,J]=dde_apply({'dde_',point1.kind,'_dot'},point1,point2,varargin{:});
end
