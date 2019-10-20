function hopf=p_tohopf(point,varargin)
%% convert point to Hopf bifurcation point
% function hopf_point=p_tohopf(point,...)
% INPUT:
%	point with stability information 
%   optional 'excludefreqs': frequency to be excluded from consideration (or
%   used otherwise, depending on input point type)
%   'funcs' problem functions
% OUTPUT:
%	hopf_point uncorrected starting guess for hopf point
%
% If point is of type 'hoho' excludefreqs can be string 'omega1' or 'omega2': then
% omega1 or omega2 will be selected as Hopf frequency, if excludefreqs is value,
% it will be excluded
%
%
% $Id: p_tohopf.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
%% for backward compatibility re-organize arguments
args=varargin;
if ~isfield(point,'x')
    funcs=point;
    point=args{1};
    args=[args(2:end),{'funcs',funcs}];
end
if isnumeric(args{1})
    args=[{'excludefreqs'},args];
    args=args(:)';
end
hopf=dde_apply({'dde_hopf_from_',point.kind,''},point,args{:});
end
