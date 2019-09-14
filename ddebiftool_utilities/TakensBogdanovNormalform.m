%% Takens-Bogdanov normal form for BT point encountered along fold or Hopf branch
%% Wrapper around locate_BT and nmfm_BT
function [nf,nflow,br_ref,indbif]=TakensBogdanovNormalform(funcs,branch,inds,varargin)
%% Input
%
% * funcs: problem functions
% * branch: fold or Hopf branch along which BT point was encountered
% * inds: array of two successive indices bracing BT point
%
%% Output
%
% * nf: point with normal form
% * nflow: if numerical finite differences are used then computation is
% done twice, once with higher order, once with lower order, this output is
% the result with lower order. Use the difference between nf and nflow to
% estimate the error
% * br_ref: refinement using extended system, including sinlge new point
% along branch
% * indbif: is the index in br_ref which is closest to BT point
%
% $Id: TakensBogdanovNormalform.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'print',0};
[options,pass_on]=dde_set_options(default,varargin,'pass_on'); 
type=branch.point(inds(1)).kind;
switch type
    case 'fold'
        tobif1=@(p)p_tofold(funcs,p);
    case 'hopf'
        tobif1=@(p)p_tohopf(funcs,p);
    otherwise        
        error('TakensBogdanovNormalform:type',...
            'TakensBogdanovNormalform: cannot detect BT along branch of type %s.',...
            type);
end
[btpoint,suc]=locate_BT(funcs,branch,inds(1));
if ~suc
    if options.print>0
        fprintf('TakensBogdanovNormalform: nonlinear system not successfully solved\n');
    end
    nf=[];
    nflow=[];
    br_ref=branch;
    indbif=[];
    return
elseif options.print>0
    fprintf('TakensBogdanovNormalform: BT point found\n');
end
btcodim1=tobif1(btpoint);
br_ref=branch;
if isfield(br_ref.point(1),'stability')
    btcodim1.stability=p_stabil(funcs,btcodim1,branch.method.stability);
end
btcodim1=dde_trim_point(btcodim1,br_ref.point(1));
br_ref.point=[br_ref.point(1:inds(1)),btcodim1,br_ref.point(inds(2):end)];
indbif=inds(2);
nmfm_compute=@(f,pt,varargin)nmfm_bt(f,pt,varargin{:});
[nf,nflow]=nmfm_wrap_findiff(funcs,btpoint,nmfm_compute,pass_on{:});
end