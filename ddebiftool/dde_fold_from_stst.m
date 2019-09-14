function fold=dde_fold_from_stst(stst,varargin)
%% convert steady state stst into approximate fold point (not corrected)
%
% optional arguments: 
%
% * 'funcs': r.h.s. functions structure to compute char. matrix if needed
%
% $Id: dde_fold_from_stst.m 308 2018-10-28 15:08:12Z jansieber $
%%
default={'funcs',[]};
options=dde_set_options(default,varargin,'pass_on');
fold=dde_fold_create('point',stst,'stability',stst.stability);
if ~isempty(stst.nvec) && strcmp(stst.flag,'fold')
    %% Fold eigenvector has been computed during normal form computations
    fold.nvec=stst.nvec;
    fold.v=fold.nvec.q;
    fold.v=fold.v/norm(fold.v);
    fold.nmfm=stst.nmfm;
    return
end
%% check if eigenvectors are provided, 
% if not, compute minimal eigenvalue of charactaristic matrix
if isempty(fold.stability) || ~isfield(fold.stability,'v')
    if isempty(options.funcs)
        error('dde_fold_from_stst:arguments',...
            'dde_fold_from_stst: eigenvectors not present in stst and r.h.s. not provided');
    end
    D=ch_matrix(options.funcs,fold.x,fold.parameter,0);
    [E1,E2]=eig(D);
    [i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
    fold.v=real(E1(:,i2));
else
   [i1,i2]=min(abs(fold.stability.l1)); %#ok<ASGLU> 
   fold.v=real(fold.stability.v(:,i2));
end
end
