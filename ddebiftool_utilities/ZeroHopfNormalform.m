%% Zero-Hopf normal form for zero-Hopf point encountered along fold or Hopf branch
%% Wrapper around nmfm_codimension2_nf
function [nf,nflow,br_ref,indbif]=ZeroHopfNormalform(funcs,branch,inds,varargin)
%% Input
%
% * funcs: problem functions
% * branch: fold or Hopf branch along which zero-Hopf point was encountered
% * inds: array of two successive indices bracing zero-Hopf point
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
% * indbif: is the index in br_ref which is closest to zero-Hopf point
%
% $Id: ZeroHopfNormalform.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'print',0};
[options,pass_on]=dde_set_options(default,varargin,'pass_on'); 
type=branch.point(inds(1)).kind;
switch type
    case 'fold'
        tobif1=@(p)p_tofold(funcs,p);
    case 'hopf'
        tobif1=@(p)p_tohopf(funcs,p,'includehopf',true);
    otherwise        
        error('ZeroHopfNormalform:type',...
            'ZeroHopfNormalform: cannot detect zero-Hopf along branch of type %s.',...
            type);
end
[zeho,suc]=locate_zeho(funcs,branch,inds(1),varargin{:});
if ~suc
    if options.print>0
        fprintf('ZeroHopfNormalform: nonlinear system not successfully solved\n');
    end
    nf=[];
    nflow=[];
    br_ref=branch;
    indbif=[];
    return
elseif options.print>0
    fprintf('ZeroHopfNormalform: zero-Hopf point found\n');
end
zehocodim1=tobif1(zeho);
br_ref=branch;
if isfield(br_ref.point(1),'stability')
    zehocodim1.stability=p_stabil(funcs,zehocodim1,branch.method.stability);
end
zehocodim1=dde_trim_point(zeho,br_ref.point(1));
br_ref.point=[br_ref.point(1:inds(1)),zehocodim1,br_ref.point(inds(2):end)];
indbif=inds(2);
nmfm_compute=@(f,pt,varargin)nmfm_zeho(f,pt,varargin{:});
[nf,nflow]=nmfm_wrap_findiff(funcs,p_tozeho(zeho),nmfm_compute,pass_on{:});
end
% type=branch.point(inds(1)).kind;
% switch type
%     case 'fold'
%         detect='hoze';
%     case 'hopf'
%         detect='zeho';
%     otherwise
%         error('ZeroHopfNormalform:type',...
%             'ZeroHopfNormalform: cannot detect Zero-Hopf along branch of type %s.',...
%             type);
% end
% nmfm_compute=@(f,pt,varargin)nmfm_zeho(f,p_tozeho(pt),varargin{:});
% [nf,nflow,br_ref,indbif]=...
%     nmfm_codimension2_nf(funcs,branch,inds,detect,nmfm_compute,varargin{:});
% type=branch.point(inds(1)).kind;
% switch type
%     case 'fold'
%         detect='hoze';
%     case 'hopf'
%         detect='zeho';
%     otherwise
%         error('ZeroHopfNormalform:type',...
%             'ZeroHopfNormalform: cannot detect Zero-Hopf along branch of type %s.',...
%             type);
% end
% nmfm_compute=@(f,pt,varargin)nmfm_zeho(f,p_tozeho(pt),varargin{:});
% [nf,nflow,br_ref,indbif]=...
%     nmfm_codimension2_nf(funcs,branch,inds,detect,nmfm_compute,varargin{:});
%end