%% Detect codimension-2 points along branch of fold points and compute all normal forms
%% Input
%
% * funcs: problem functions
% * branch: fold branch
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
% * testfuncs: structure containing fields 'cusp', 'zeho', 'BT' with
% sign-changing test functions
%
% If an entry in bifpoints is empty, then detection or coputation have failed.
%%
function [bifpoints,indices,branch,testfuncs]=FoldCodimension2(funcs,branch,varargin)
%
% $Id: FoldCodimension2.m 309 2018-10-28 19:02:42Z jansieber $
%
%% 
default={'max_nferror',1e-3,'print',1};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch=replace_branch_pars(branch,[],pass_on);
%% determine stability and normal form coefficients
if options.print>0
    fprintf('FoldCodimension2: calculate stability if not yet present\n');
end
[nunst,dom,triv,branch.point]=GetStability(branch,'funcs',funcs,...
    'exclude_trivial',true,varargin{:}); %#ok<ASGLU>
if options.print>0
    fprintf('FoldCodimension2: calculate fold normal form coefficients\n');
end
[b,blow,branch]=NormalformCoefficients(funcs,branch,'print',options.print-1);
%% check for large discrepancies (if finite differences are used)
b_error_large=find(abs(b-blow)./max(abs(b),1)>options.max_nferror);
if ~isempty(b_error_large)
    warning('FoldComdimension2:b',...
        ['FoldCodimension2: large error for normal form coefficient in points\n',...
        '%s.'],num2str(b_error_large));
end
%% Find all intervals for special points
bchange=diff(sign(b))~=0;                    % for cusp
nby1=abs(diff(nunst(:)'))==1;                % for Takens-Bogdanov
nby2=abs(diff(nunst(:)'))==2;                % for zero-Hopf
nbymore=abs(diff(nunst(:)'))>2;
bchange=bchange & (~nby1) & (~nbymore);  % fold coefficient is singular at BTs
detect=cell(1,length(branch.point));
rep=@(f,sel)repmat({f},1,sum(sel));
detect(bchange)=rep(@CuspNormalform,bchange);
detect(nby1)=rep(@TakensBogdanovNormalform,nby1);
detect(nby2)=rep(@ZeroHopfNormalform,nby2);
detect(nbymore)=rep(@nmfm_dummy,nbymore);
if options.print>0
    fprintf('FoldCodimension2: (provisional) ');
    nmfm_printbif_type('cusp',bchange);
    nmfm_printbif_type('Takens-Bogdanov',nby1);
    nmfm_printbif_type('Zero-Hopf',nby2);
    nmfm_printbif_type('unknown',nbymore);
    fprintf(' detected.\n');
end    

%% Detect bifurcation points and compute their normal form
[bifpoints,biflow,branch,indices]=br_insert(funcs,branch,detect,pass_on{:},'print',options.print);
%% check if for any bifurcations normal form coefficients seem wrong
s2a=@(s)structfun(@(x)x,s.nmfm);
nfdiffs=cellfun(@(p,pl)norm(s2a(p)-s2a(pl),'inf'),bifpoints,biflow);
nf_error_large=find(nfdiffs>options.max_nferror);
if ~isempty(nf_error_large)
    warning('FoldComdimension2:b',...
        ['FoldCodimension2: large error for normal form coefficient in codim2 points\n',...
        '%s:\n%s'],num2str(nf_error_large),num2str(nfdiffs(nf_error_large)));
end
%% (Re)compute test functions if required
if nargout>3
    [b,blow]=NormalformCoefficients(funcs,branch);
    testfuncs.cusp=[b;blow];
    testfuncs.bt=phi_f_bt(funcs,branch.point);
    imagthresh=branch.method.bifurcation.imagthreshold; % threshold for treating a root as complex
    isimag=@(x)x>imagthresh;
    testfuncs.zeho=arrayfun(...
        @(p)nmfm_smrp(funcs, p,branch.method.stability,'remove_omega',false,'threshold',isimag),...
        branch.point);
end
end
