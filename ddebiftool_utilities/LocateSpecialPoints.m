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
%% Simple wrapper around branch specific functions
types=struct(...
    'stst',@StstCodimension1,...
    'fold',@FoldCodimension2,...
    'hopf',@HopfCodimension2);
btype=branch.point(1).kind;
if ~isfield(types,btype)
    warning('LocateSpecialPoints:type',...
        'LocateSpecialPoints for type %s not implemented', btype);
    testfuncs=[];
    indices=[];
    biftype={};
    return
end
routine=types.(btype);
[bifpoints,indices,branch,testfuncs]=routine(funcs,branch,varargin{:});
%% insert bifpoints into refined branch
%#ok<*SFLD>
branch=br_flag(branch);biftype=cell(1,length(indices));
for i=1:length(indices)
    biftype{i}=bifpoints{i}.kind;
    branch.point(indices(i)).flag=bifpoints{i}.kind;
    branch.point(indices(i)).nmfm=structmerge(branch.point(indices(i)).nmfm,bifpoints{i}.nmfm);
    % below: only copy codim2 normal form nvec if nvec is empty. codim2
    % don't need previous eigenvectors while codim1 normal forms need it
    % for continuity of test function. However, Hopf-Hopf detection of
    % trivial eigenvalues profits from nvec (contains omega1,omega2)
    if isfield(bifpoints{i},'nvec') && isempty(branch.point(indices(i)).nvec)
        branch.point(indices(i)).nvec=bifpoints{i}.nvec;
    end
%     if isfield(bifpoints{i},'omega')
%         branch.point(indices(i)).omega=bifpoints{i}.omega;
%     end
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
