%% codimension-2 normal form for special point encountered along codimension-1 branch
%%
function [nf,nflow,br_ref,indbif]=...
    nmfm_codimension2_nf(funcs,branch,inds,detect,nmfm_compute,varargin)
%% Input
%
% * funcs: problem functions
% * branch: codim1 branch along which codim-2 poitn was encountered
% * inds: array of two successive indices bracing codimension-2 point
% * detect: monitoring function of type res=detect(p) for points of branch,
% or string ('hoho', 'zeho','hoze')
% * nmfm_compute: function of type newpoint=nmfm_type(funcs,point) where
% type is hoho, zeho, BT, cusp etc
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
% * indbif: is the index in br_ref which is closest to the sign change of
% the monitoring function detect.
%
% $Id: nmfm_codimension2_nf.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'print',0};
[options,pass_on]=dde_set_options(default,varargin,'pass_on'); 
[br_ref,indbif,indmap,notcorrected]=...
    br_bisection(funcs,branch,inds,detect,pass_on{:},'print',options.print-1); %#ok<ASGLU>
if options.print>0
    if ~isempty(notcorrected) && notcorrected(indbif)
        fprintf(['nmfm_codinmension2: ',...
            'Newton correction failed in proposed bifurcation point,\n',...
            'using secant approximation\n']);
    end
end
if isfinite(indbif) && indbif>=1 && indbif<=length(br_ref.point)
    pt=br_ref.point(indbif);
    [nf,nflow]=nmfm_wrap_findiff(funcs,pt,nmfm_compute,pass_on{:});
else
    nf=br_ref.point([]);
    nflow=nf;
end
end

