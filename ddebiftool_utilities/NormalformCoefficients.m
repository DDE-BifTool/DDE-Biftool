%% Normal form coefficient a in fold point along fold or Hopf branch
%% Input
%
% * funcs: problem functions
% * branch: fold branch or array of points of kind 'fold' or 'hopf'
% * optional name-value pair: 
% * 'use_nullvectors' (boolean, default false, for hopf only) use nullvectors from
% previous point along branch to compute new nullvectors by bordering
% (otherwise, svd is used)
%
%% Output
%
% * coeff: array of normal form coefficients along branch
% * coefflow: equals b if funcs.sys_mfderi is provided, otherwise (for
% finite-difference approximation) a lower-order estimate of coeff
% coefficients. Use coeff-coefflow to estimate the accuracy of coeff. If this is
% large, finite-differences are problematic. If this is small it may(!) be
% ok.
% * branch: branch where points have nmfm field attached
% * updated: logical array indicating which normal forms have been
% recomputed (if nmfm is present and non-empty and optinoal input recompute
% is false, only non-present normal forms get recomputed.)
%%
function [coeff,coefflow,branch,updated]=NormalformCoefficients(funcs,branch,varargin)
%%
%
% $Id: NormalformCoefficients.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'print',0,'recompute',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isfield(branch,'point')
    pt=branch.point;
else
    pt=branch;
end
switch pt(1).kind
    case 'hopf'
        use_nullvectors=false;
        nffunc=@(fun,point,np,varargin)nmfm_hopf(fun,point,'nullpoint',np,varargin{:});
        cfield='L1';
    case 'fold'
        use_nullvectors=true;
        nffunc=@(fun,point,np,varargin)nmfm_fold(fun,point,'nullpoint',np,varargin{:});
        cfield='b';
end
npt=length(pt);
nullpoint=[];
coeff=NaN(1,npt);
coefflow=coeff;
updated=false(1,npt);
for i=npt:-1:1
    currpt=pt(i);
    if ~options.recompute && isfield(currpt,'nmfm') && ~isempty(currpt.nmfm)
        nf=currpt;
        nflow=currpt;
    else
        nmfm_compute=@(fun,point,varargin)nffunc(fun,point,nullpoint,varargin{:});
        [nf,nflow]=nmfm_wrap_findiff(funcs,currpt,nmfm_compute,pass_on{:});
        updated(i)=true;
    end
    if isfield(nf.nmfm,cfield)
        coeff(i)=nf.nmfm.(cfield);
    else
        coeff(i)=NaN;
    end
    if isfield(nflow.nmfm,cfield)
        coefflow(i)=nflow.nmfm.(cfield);
    else
        coefflow(i)=NaN;
    end
    if options.print>0
        fprintf('Normalform Coefficients: pt %d of %d on %s branch: %s=%g\n',...
            i,npt,pt(1).kind,cfield,coeff(i));
    end
    if use_nullvectors && isfield(nf.nvec,'q') && all(isfinite(nf.nvec.q)) && all(isfinite(nf.nvec.q)) 
        nullpoint=nf;
    end
    ptnew(i)=nf;
end
if isfield(branch,'point')
    branch.point=ptnew;
else
    branch=ptnew;
end
end
