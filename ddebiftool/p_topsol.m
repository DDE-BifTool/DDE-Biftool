function [psol,stpcond]=p_topsol(point,varargin)
%% create starting guess for periodic solution derived from point
% function [psol_point,stpcond]=p_topsol(funcs,point,ampl,col_degree,nr_int)
% INPUT:
%	point: point of kind hopf, psol or hcli (initially)
%
% optional arguments
% 
%	'radius' amplitude of periodic solution guess
%       % from Hopf:
%	'degree' piecewise polynomial degree for periodic solution
%	'intervals' number of intervals for mesh
%       % for period doubling:
%	this information, but it is needed for branching off
%
% OUTPUT:
%	psol (uncorrected) starting guess for periodic solution derived from point
%	stpcond steplength condition for use in correction of guess

% (c) DDE-BIFTOOL v. 2.02, 30/06/2002
%
% $Id: p_topsol.m 308 2018-10-28 15:08:12Z jansieber $
%% for backward compatibility re-organize arguments
args=varargin;
if ~isfield(point,'kind')
    funcs=point;
    point=args{1};
    args=[args(2:end),{'funcs',funcs}];
end
if isnumeric(args{1})
    argsfixed=[{'radius','degree','intervals'};args(1:3)];
    args=[argsfixed(:)',args(4:end)];
end
%% call source specific function (hopf, psol and hcli supported initially)
[psol,stpcond]=dde_apply({'dde_psol_from_',point.kind,''},point,args{:});
end
