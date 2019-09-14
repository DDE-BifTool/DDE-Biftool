%% Detect bifurcation points along branch of steady states and compute all normal forms
%% Input
%
% * funcs: problem functions
% * branch: stst branch
% * optional name-value pairs:
% * 'max_nferror' (default=|1e-3|) warn if difference between different
% approximation orders is larger than this (useful for finite-difference
% approximation only)
% * all other name-value pairs are passed on as fields for branch fields
%% Output
% * bifpoints: 1 x ns cell array of special points
% * indices: 1 x ns array of pointers into branch.point at which special
% points are inserted (as ordinary fold points)
% * branch: updated branch with inserted special points
% * testfuncs: structure containing fields 'fold', 'hopf' with
% sign-changing test functions
%
% If an entry in bifpoints is empty, then detection or computation have failed.
%%
function [bifpoints,indices,branch,testfuncs]=StstCodimension1(funcs,branch,varargin)
%
% $Id: StstCodimension1.m 309 2018-10-28 19:02:42Z jansieber $
%
%% 
default={'max_nferror',1e-3,'print',1};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch=replace_branch_pars(branch,[],pass_on);
%% determine stability, normal form coefficients and det(Delta(lambda))
foldfunc=@(p)det(ch_matrix(funcs,p.x,p.parameter,0));
foldbr=@(br)arrayfun(foldfunc,br.point);
if options.print>0
    fprintf('StstCodimension1: calculate stability if not yet present\n');
end
[nunst,dom,triv,branch.point]=GetStability(branch,'funcs',funcs); %#ok<ASGLU>
fold=foldbr(branch);
%% Find all intervals for special points
nby1=abs(diff(sign(fold(:)')))~=0;        % for fold
nby2=abs(diff(nunst(:)'))==2;             % for Hopf
nbymore=abs(diff(nunst(:)'))>2;
detect=cell(1,length(branch.point));
rep=@(f,sel)repmat({f},1,sum(sel));
detect(nby1)=rep(@(funcs,branch,inds,varargin)...
    Codim1Normalform(funcs,branch,inds,'fold',varargin{:}),nby1);
detect(nby2)=rep(@(funcs,branch,inds,varargin)...
    Codim1Normalform(funcs,branch,inds,'hopf',varargin{:}),nby2);
detect(nbymore)=rep(@nmfm_dummy,nbymore);
if options.print>0
    fprintf('StstCodimension1: (provisional) ');
    nmfm_printbif_type('fold',nby1);
    nmfm_printbif_type('Hopf',nby2);
    nmfm_printbif_type('unknown',nbymore);
    fprintf(' detected.\n');
end    
%% Detect bifurcation points and compute their normal form
[bifpoints,biflow,branch,indices]=br_insert(funcs,branch,detect,pass_on{:},...
    'print',options.print);
%% check if for any bifurcations normal form coefficients seem wrong
s2a=@(s)structfun(@(x)x,s.nmfm);
nfdiffs=cellfun(@(p,pl)norm(s2a(p)-s2a(pl),'inf'),bifpoints,biflow);
nf_error_large=find(nfdiffs>options.max_nferror);
if ~isempty(nf_error_large)
    warning('StstComdimension1:L1etc',...
        ['StstCodimension1: large error for normal form coefficient(s) in bifurcation points\n',...
        '%s:\n%s'],num2str(nf_error_large),num2str(nfdiffs(nf_error_large)));
end
%% (Re)compute test functions if required
if nargout>3
    testfuncs.fold=foldbr(branch);
    imagthresh=branch.method.bifurcation.imagthreshold; % threshold for treating a root as complex
    isimag=@(x)x>imagthresh;
    testfuncs.hopf=arrayfun(...
        @(p)nmfm_smrp(funcs, p,branch.method.stability,'threshold',isimag),...
        branch.point);
end
end
