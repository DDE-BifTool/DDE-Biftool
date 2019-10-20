function hopf=dde_hopf_from_stst(stst,varargin)
%% convert steady state stst into approximate hopf point (not corrected)
%
% optional arguments: 
%
% * 'excludefreqs': list of frequencies to exclude
% * 'funcs': r.h.s. functions structure to compute eigenvectors if needed
% * 'method': parameters for stability computation
% $Id: dde_hopf_from_stst.m 315 2019-01-29 19:42:21Z jansieber $
%%
default={'funcs',[],'excludefreqs',[],'includehopf',false,...
    'method',getfield(df_mthod('stst'),'stability')};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
hopf=setfield(dde_hopf_create(stst),'stability',stst.stability);
excludefreqs=abs(options.excludefreqs(:).');
if ~isempty(hopf.omega) && ~options.includehopf
    excludefreqs=[abs(hopf.omega),excludefreqs];
elseif isfield(stst,'nvec') && isfield(stst,'flag') && ~isempty(stst.nvec) && strcmp(stst.flag,'hopf')
    %% Hopf eigenvector and value has been computed during normal form computations
    hopf.nvec=stst.nvec;
    hopf.omega=abs(hopf.nvec.omega);
    hopf.v=hopf.nvec.q;
    hopf.v=hopf.v/norm(hopf.v);
    hopf.nmfm=stst.nmfm;
    return
end
%% check if eigenvectors are provided, if not, compute
if isempty(hopf.stability) || ~isfield(hopf.stability,'v')
    if isempty(options.funcs)
        error('dde_hopf_from_stst:arguments',...
            'dde_hopf_from_stst: eigenvectors not present in stst and r.h.s. not provided');
    end
    hopf.stability=p_stabil(options.funcs,hopf,options.method,pass_on{:});
end
evals=hopf.stability.l1;
selimp=find(imag(evals)>0);
%% remove frequencies to be excluded
if ~isempty(excludefreqs)
    ind=dde_match_complex(excludefreqs(:),imag(evals(selimp)));
    selimp(ind)=[];
end
if isempty(selimp)
    error('dde_hopf_from_stst:eigenvalues',['dde_hopf_from_stst: ',...
        'no good pair of complex roots found.']);
end
%% select eigenvalue closest to imaginary axis
evals=evals(selimp);
[i1,i2]=min(abs(real(evals))); %#ok<ASGLU>
hopf.omega=imag(evals(i2));
%% find corresponding eigenvector
if isfield(hopf.stability,'v')
    hopf.v=hopf.stability.v(:,selimp(i2));
elseif ~isempty(options.funcs)
    D=ch_matrix(options.funcs,hopf.x,hopf.parameter,1i*hopf.omega);
    [E1,E2]=eig(D);
    [i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
    hopf.v=E1(:,i2);
else
    error('dde_hopf_from_stst:arguments',...
        'dde_hopf_from_stst: eigenvectors not present in stst and r.h.s. not provided');
end
end
