%% Detect special points along branch of points and compute all normal forms
%% Input
%
% * funcs: problem functions
% * branch: branch of solutions (at the moment supported: hopf, fold, stst)
% * optional name-value pairs:
% * 'max_nferror' (default=|1e-3|) warn if difference between different
% approximation orders is larger than this (useful for finite-difference
% approximation only)
% * all other name-value pairs are passed on as fields for branch fields
%% Output
% * branch: updated branch with inserted special points
% * testfuncs: structure containing fields with branch specific test functions
% * indices: locations of special points in branch.point
% * biftype: cell array with names of bifurcations
%
% If an entry in bifpoints is empty detection or computation have failed.
%%
function [branch,testfuncs,indices,biftype]=LocateSpecialPoints(funcs,branch,varargin)
%
% $Id: LocateSpecialPoints.m 309 2018-10-28 19:02:42Z jansieber $
%
default={'pointtype_list',pointtype_list()};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% Simple wrapper around branch specific functions
testfuncs=[];
indices=[];
biftype={};
if isempty(branch.point)
    return
end
if isfield(funcs,'kind')
    kind=funcs.kind;
else
    kind=branch.point(1).kind;
end
Kind=[upper(kind(1)),kind(2:end)];
codim=options.pointtype_list.codim.(kind);
routine=[Kind,'Codimension',num2str(codim+1)];
if ~exist(routine,'file')
    warning('LocateSpecialPoints:type',...
        'LocateSpecialPoints for type %s not implemented', kind);
    return
end
[bifpoints,indices,branch,testfuncs]=feval(routine,funcs,branch,pass_on{:});
%% insert bifpoints into refined branch
%#ok<*SFLD>
biftype=cell(1,length(indices));
for i=1:length(indices)
    biftype{i}=bifpoints{i}.kind;
    if isfield(branch.point(indices(i)),'flag')
        branch.point(indices(i)).flag=bifpoints{i}.kind;
    end
    if isfield(branch.point(indices(i)),'nmfm')
        branch.point(indices(i)).nmfm=structmerge(branch.point(indices(i)).nmfm,bifpoints{i}.nmfm);
    end
    % below: only copy codim2 normal form nvec if nvec is empty. codim2
    % don't need previous eigenvectors while codim1 normal forms need it
    % for continuity of test function. However, Hopf-Hopf detection of
    % trivial eigenvalues profits from nvec (contains omega1,omega2)
    if isfield(bifpoints{i},'nvec') && isempty(branch.point(indices(i)).nvec)
        branch.point(indices(i)).nvec=bifpoints{i}.nvec;
    end
end
end
function s=structmerge(varargin)
s=struct();
for k=1:length(varargin)
    if isempty(varargin{k})
        continue
    end
    names=fieldnames(varargin{k});
    for i=1:length(names)
        s.(names{i})=varargin{k}.(names{i});
    end
end
end
