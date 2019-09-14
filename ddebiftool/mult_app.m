function [mu,eigenfuncs]=mult_app(funcs,period,profile,mesh,degree,rho,max_number,col,par,d_ac)
%% find Floquet multipliers & modes of periodic orbit
% (wrapper around dde_psol_eig kept for backward compatibility)
% function [mu,eigenfuncs]=mult_app(period,profile,mesh,degree,rho,max_number,col,par,d_ac)
% INPUT: 
%   funcs problem functions
%	period period of solution
%	profile periodic solution profile
%       mesh periodic solution mesh (if empty, mesh is assumed uniform)
%	degree degree of piecewise polynomials
%       rho keep multipliers with modulus >= rho
%	max_number keep at most max_number multipliers 
%	col collocation points
%       par current parameter values in R^p
%       d_ac (only for state-dependent delays) tau<d_ac is treated as 
%             tau<0 (stability is not computed)
% OUTPUT:
%       mu approximations of requested multipliers
%       eigenfuncs (if requested) corresponding modes

%  (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: mult_app.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
if nargin<10
    d_ac=1e-8;
end
%% if mesh is empty assume equidistant mesh:
if isempty(mesh)
    mesh=0:1/(size(profile,2)-1):1;
end
if nargout>1
    geteigenfuncs=true;
else
    geteigenfuncs=false;
end
psol=dde_psol_create('mesh',mesh,'profile',profile,'period',period,'degree',degree,'parameter',par);
stab=dde_psol_eig(funcs,psol,'collocation_parameters',col,'delay_accuracy',d_ac,'minimal_modulus',rho,...
    'max_number_of_eigenvalues',max_number,'geteigenfuncs',geteigenfuncs);
mu=stab.mu;
if geteigenfuncs
    eigenfuncs=reshape([stab.eigenfuncs.profile],numel(profile),length(stab.eigenfuncs));
end
end
