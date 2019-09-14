%% Number of unstable eigenvalues/Floquet multipliers along branch
%%
function [nunst,dom,triv_defect,points]=GetStability(branch,varargin)
%% Inputs
%
% * |branch|: branch along which number of unstable eigenvalues/Floquet
% multipliers is determined
%
%  Important name-value pair inputs
% 
% * |'exclude_trivial'| (logical): exclude trivial eigenvalues (what they are
%                            depends on the type of points)
% * |'locate_trivial'| (function): e.g. @(p)[1;-1] for removing 1 and -1
%                            for period doubling orbits (overwrites standard location)
% * |'funcs'|:   set of functions for computing branch, needed if stability 
% information is not yet present and needs to be calculated
% * |'method'|: which method is used for comutation (if not present, a
% standard method is generated with |df_mthod|)
% * |'recompute'| (default |false|):  force recomputation of stability even
% if it is present.
%
%% Outputs
% * |nunst| (vector of integers): number of unstable eigenvalues
% * |dom| (vector of real/complex): dominant eigenvalue (closest unstable to
%                       bifurcation if exists, otherwise, closest stable)
% * |triv_defect| (vector of real): distance of eigenvalue approximating
%                               trivial eigenvalue from its true value
% * |points| (struct array): array of points (same as in input branch but
% with stability info if not present on input)
%
% $Id: GetStability.m 369 2019-08-27 00:07:02Z jansieber $
%
%% process options
defaults={'exclude_trivial',false,'locate_trivial',[],'funcs',[],'critical',false,...
    'points',[],'method',[],'recompute',false,...
    'pointtype_list',@pointtype_list};
[options,pass_on]=dde_set_options(defaults,varargin,'pass_on');
%% check if stability has been computed already & compute if necessary
if isfield(branch,'point')
    points=branch.point;
    if ~isempty(options.points) % select only some points if desired
        points=branch.point(options.points);
    end
else
    points=branch;
end
%% for periodic orbit bifurcations extract solution component
kind=points(1).kind;
funcs=options.funcs;
if ~isempty(options.funcs) && isfield(options.funcs,'kind') && ...
    isfield(options.funcs,'get_comp') && isfield(options.funcs,'userfuncs')
    kind=options.funcs.kind;
    fullpoints=points;
    points=options.funcs.get_comp(points,'solution_for_stability');
    funcs=options.funcs.userfuncs;
else
    fullpoints=[];
end
%% detrmine stability
if options.recompute || ~isfield(points,'stability') ||...
        any(arrayfun(@(x)isempty(x.stability),points))
    %% create method for calculating stability
    if ~isempty(options.method)
        mth=options.method;
    elseif isfield(branch,'method') && isfield(branch.method,'stability')
        mth=branch.method;
        mth.stability=dde_set_options(mth.stability,pass_on,'pass_on');
    else
        mth=df_mthod(funcs,points(1).kind);
        if isfield(mth,'stability')
            mth.stability=dde_set_options(mth.stability,pass_on,'pass_on');
        end
    end
    for i=1:length(points)
        if isfield(mth,'stability') && ...
                (options.recompute || ~isfield(points(i),'stability') || isempty(points(i).stability))
            stab=p_stabil(funcs,points(i),mth.stability,pass_on{:});
            points(i).stability=stab;
        end
    end
end
%% Count number of unstable eigenvalues
% (removing trivial ones as indicated by |locate_trivial| input or standard
% theory)
np=length(points);
nunst=NaN(np,1);
dom=nunst;
triv_defect=zeros(np,1);
% get definitions of stability and trivial eigenvalues
evfcn=options.pointtype_list(pass_on{:}); 
if ~isempty(options.locate_trivial)
    % overwrite location of trivial eigenvalues
    evfcn.triv.(kind)=options.locate_trivial;
end
for i=1:np
    p=points(i);
    %% check if point is special point
    if isfield(p,'flag') && ~isempty(p.flag)
        type=p.flag;
    else
        type=kind;
    end
    %% exclude trivial eigenvalues
    excl=evfcn.triv.(type)(p);
    ev=evfcn.getev.(type)(p);
    if options.exclude_trivial
        for k=1:length(excl)
            [dist,ind]=min(abs(ev-excl(k)));
            ev(ind)=[];
            if ~isempty(dist)
                triv_defect(i)=max(dist,triv_defect(i));
            end
        end
    end
    %% count unstable eigenvalues
    unstsel=evfcn.stab.(type)(ev)>=0;
    nunstloc=sum(unstsel);
    if isempty(nunstloc)
        nunst(i)=0;
    else
        nunst(i)=nunstloc;
    end
    %% check which other eigenvalue is dominant (critical means largest)
    if nunst(i)>0 && ~options.critical
        [dum,ind]=min(evfcn.stab.(type)(ev(unstsel))); %#ok<ASGLU>
        dom(i)=ev(ind);
        if imag(ev(ind))<0
            dom(i)=conj(dom(i));
        end
    else
        [dum,ind]=min(abs(evfcn.stab.(type)(ev)));
        if ~isempty(dum)
            dom(i)=ev(ind);
            if imag(ev(ind))<0
                dom(i)=conj(dom(i));
            end
        end
    end
end
%% for psol bifurcations, remove extra fields and re-insert extended components
if ~isempty(fullpoints)
    points=arrayfun(@(x,y)setfield(x,'stability',y.stability),fullpoints,points); %#ok<SFLD>
end
end