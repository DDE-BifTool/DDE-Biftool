%% Detect codimension-2 points along branch of Hopf points and compute all normal forms
%% Input
%
% * funcs: problem functions
% * branch: Hopf branch
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
% If an entry in bifpoints is empty detection or computation have failed.
%%
function [bifpoints,indices,branch,testfuncs]=HopfCodimension2(funcs,branch,varargin)
%
% $Id: HopfCodimension2.m 309 2018-10-28 19:02:42Z jansieber $
%
%% 
default={'max_nferror',1e-3,'print',1,'stabilityfield','l0'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch=replace_branch_pars(branch,[],pass_on);
%% determine stability, normal form coefficients and det(Delta(lambda))
zehofunc=@(p)det(ch_matrix(funcs,p.x,p.parameter,0));
zehobr=@(br)arrayfun(zehofunc,br.point);
if options.print>0
    fprintf('HopfCodimension2: calculate stability if not yet present\n');
end
[nunst,dom,triv,branch.point]=GetStability(branch,'funcs',funcs,...
    'exclude_trivial',true,'stabilityfield',options.stabilityfield); %#ok<ASGLU>
zeho=zehobr(branch);
if options.print>0
    fprintf('HopfCodimension2: calculate L1 coefficients\n');
end
[L1,L1low,branch]=NormalformCoefficients(funcs,branch,pass_on{:},'print',options.print-1);
%% check for large relative discrepancies (if finite differences are used)
b_error_large=find(abs(L1-L1low)./max(abs(L1),1)>options.max_nferror);
if ~isempty(b_error_large)
    warning('HopfComdimension2:L1',...
        ['HopfCodimension2: large error for L1 in points\n',...
        '%s.'],num2str(b_error_large));
end
%% Find all intervals for special points
L1change=diff(sign(L1))~=0;               % for generalized Hopf
nby1=abs(diff(sign(zeho(:)')))~=0;        % for zero-Hopf
nby2=abs(diff(nunst(:)'))==2;             % for Hopf-Hopf
omega=arrayfun(@(x)x.omega,branch.point);
omsign=diff(sign(omega))~=0;              % for Takens-Bogdanov
omzero=omega(1:end-1)==0|omega(2:end)==0;
nbymore=abs(diff(nunst(:)'))>2;
% L1 coefficient is singular at BTs and zehos
L1change=L1change & (~omsign) & (~omzero) & (~nbymore) & (~nby1); 
detect=cell(1,length(nunst));
rep=@(f,sel)repmat({f},1,sum(sel));
detect(L1change)=rep(@GeneralizedHopfNormalform,L1change);
detect(nby1)=rep(@ZeroHopfNormalform,nby1);
detect(nby2)=rep(@HopfHopfNormalform,nby2);
detect(omsign)=rep(@TakensBogdanovNormalform,omsign);
detect(nbymore)=rep(@nmfm_dummy,nbymore);
if options.print>0
    fprintf('HopfCodimension2: (provisional) ');
    nmfm_printbif_type('gen. Hopf',L1change);
    nmfm_printbif_type('Takens-Bogdanov',omsign);
    nmfm_printbif_type('Zero-Hopf',nby1);
    nmfm_printbif_type('Hopf-Hopf',nby2);
    nmfm_printbif_type('unknown',nbymore);
    fprintf(' detected.\n');
end    
%% Detect bifurcation points and compute their normal form
[bifpoints,biflow,branch,indices]=br_insert(funcs,branch,detect,pass_on{:},...
    'print',options.print,'stabilityfield',options.stabilityfield);
%% check if for any bifurcations normal form coefficients seem wrong
% problem when nmfm field contains vectors, needs to be fixed 
% s2a=@(s)structfun(@(x)x,s.nmfm);
% nfdiffs=cellfun(@(p,pl)norm(s2a(p)-s2a(pl),'inf'),bifpoints,biflow);
% nf_error_large=find(nfdiffs>options.max_nferror);
% if ~isempty(nf_error_large)
%     warning('HopfComdimension2:L2etc',...
%         ['HopfCodimension2: large error for normal form coefficient(s) in codim2 points\n',...
%         '%s:\n%s'],num2str(nf_error_large),num2str(nfdiffs(nf_error_large)));
% end
%% (Re)compute test functions if required
if nargout>3
    if sum(omsign)>0
        warning('off','NMFM_HOPF:omega'); % at Takens-Bogdanov L1 is singular
    end
    L1old=L1;
    L1lowold=L1low;
    [L1,L1low,branch,updated]=NormalformCoefficients(funcs,branch,pass_on{:},'print',options.print-1);
    L1(~updated)=L1old;
    L1low(~updated)=L1lowold;
    if sum(omsign)>0
        warning('on','NMFM_HOPF:omega');
    end
    testfuncs.genh=[L1;L1low];
    testfuncs.bt=arrayfun(@(x)x.omega,branch.point);
    testfuncs.zeho=zehobr(branch);
    imagthresh=branch.method.bifurcation.imagthreshold; % threshold for treating a root as complex
    isimag=@(x)x>imagthresh;
    testfuncs.hoho=arrayfun(...
        @(p)nmfm_smrp(funcs, p,branch.method.stability,'remove_omega',true,'threshold',isimag),...
        branch.point);
end
end
