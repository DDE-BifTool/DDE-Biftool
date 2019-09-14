%% construct eigenvalue problem for Floquet multipliers using sparse matrices
% for case with only one-sided delays. This is based on eigs and permits
% optional argument 'closest' for finding Floquet multipliers close to a
% particular value.
%
% $Id: dde_psol_mult_app_sparse.m 352 2019-06-19 16:03:03Z jansieber $
%%
function [s,ef]=dde_psol_mult_app_sparse(Margs,varargin)
default={'geteigenfuncs',false,...
    'closest',[],'max_number_of_eigenvalues',20,'method',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on','method');
if isnumeric(Margs{1})
    dim=size(Margs{1},1);
else
    dim=Margs{2};
end
n_ev=min(options.max_number_of_eigenvalues,dim);
closest=num2cell(options.closest);
if isempty(closest)
    closest={'lm'};
end
if ~options.geteigenfuncs
    s=eigs(Margs{:},n_ev,closest{:});
    ef=[];
else
    [ef,s]=eigs(Margs{:},n_ev,closest{:});
    s=diag(s);
end
end
