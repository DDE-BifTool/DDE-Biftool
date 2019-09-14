%% Hopf-Hopf normal form for Hopf-Hopf point encountered along Hopf branch
%% Wrapper around nmfm_codimension2_nf
%% Input
%
% * funcs: problem functions
% * branch: Hopf branch along which Hopf-Hopf point was encountered
% * inds: array of two successive point indices bracing Hopf-Hopf point
%
%% Output
%
% * nf: point with normal form
% * nflow: if numerical finite differences are used then computation is
% done twice, once with higher order, once with lower order, this output is
% the result with lower order. Use the difference between nf and nflow to
% estimate the error
% * br_ref: bisection refinements are performed along the branch before
% normal form calculation, br_ref is branch with additional points between
% indices inds
% * indbif: is the index in br_ref which is closest to Hopf-Hopf point
%
% $Id: HopfHopfNormalform.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
function [nf,nflow,br_ref,indbif]=HopfHopfNormalform(funcs,branch,inds,varargin)
if ~strcmp(branch.point(inds(1)).kind,'hopf')
    error('HopfHopfNormalform:type',...
            'HopfHopfNormalform: cannot detect Hopf-Hopf along branch of type %s.',...
            type);
end
nmfm_compute=@(f,pt,varargin)nmfm_hoho(f,p_tohoho(pt),varargin{:});
[nf,nflow,br_ref,indbif]=...
    nmfm_codimension2_nf(funcs,branch,inds,'hoho',nmfm_compute,varargin{:});
end
