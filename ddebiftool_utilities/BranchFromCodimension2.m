%% Branch off from special points along branch
%
% This helper function creates all available information to branch off from
% codim 2 points along a codimension 1 branch (at the moment |branch| needs
% to be a branch of bifurcations of equilibria: type |fold| or |hopf|).
%
%% Inputs
%
% * |funcs|: problems description structure ,created, e.g., by |set_funcs|
% * |branch|: branch structure, including normal form information fields
% |'nvec'|, |'nmfm'| and |'flag'| in its |'point'| fields. Valid kinds of
% |'point'| are fold and hopf.
%
%% Optional name-value pairs
% If no optional arguments are given, all branching information for all
% codimension two points along the branch is computed.
%
% * |'indices'| (empty): only determine branching information for
% specific points along the branch (numbers given here)
% * |'name'| (empty): only determine branches of this type (e.g.,
% 'TorusBifurcation'). If empty, all branches are taken as found by
% |nmfm_branch_list(..)|.
% * |'selection'| (|'all'|, needed if |'name'| is non-empty). Select only
% some of the branches. Say, if arguments are
% |...,'name','hopf','selection',1)| then only the first Hopf bifurcation
% branching off is reported. In case of Hopf-Hopf points this skips the
% current branch.
% * |'print'| (2): print level. Level 2 prints out which branches have been
% discovered and convergence success information for the initial points.
%
% Other name-value pairs are passed on to called functions. Useful options
% for branch initialization used in |C1branch_from_C2point|:
%
% * |'step'| (0.01): initial step size along branch in secant distance
% * |'correc'| (true): flag controlling whether to perform Newton iteration
% * |'degree'| (3), |'intervals'| (20): discretization mesh size for
% periodic orbit bifurcation points.
%% Output
%
% Array of structures with fields
%
% * |'C1type'|: bifurcation branching off (e.g., |'TorusBifurcation'|)
% * |'C2type'|: flag of codim 2 point
% * |'funcs'|: problem structure to be used for subsequent |br_contn| call
% * |'reverse'|: hint for subsequent |br_contn| call that branch can be
% reversed
% * |'success'|: flag, true if all initial corrections were successful.
% Continuation can still be started if |success| is |false|, but may fail.
% 
% $Id: BranchFromCodimension2.m 309 2018-10-28 19:02:42Z jansieber $
%%
function C1info=BranchFromCodimension2(funcs,branch,varargin)
%%
default={'name',[],'selection','all','reverse',false,'indices',[],'print',2};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
C1info=[];
if isempty(options.indices)
    options.indices=find(~strcmp('',{branch.point.flag}));
end
for ind=1:length(options.indices)
    C2point=branch.point(options.indices(ind));
    C2type=C2point.flag;
    if isempty(options.name)
        namelist=nmfm_branch_list(C2type);
    else
        namelist=struct(...
            'name',options.name,...
            'selection',options.selection,...
            'reverse',options.reverse);
    end
    if isempty(namelist)
        if options.print>0
            fprintf(' Point %d, flag=%s: no branching instructions found.\n',...
                ind, C2type);
        end
        continue
    end
    for i=length(namelist):-1:1
        entry=namelist(i);
        if options.print>0
            fprintf('from %s at',C2type);
            for k=1:length(branch.parameter.free)
                p_ind=branch.parameter.free(k);
                fprintf(' par(%d)=%g,',p_ind,C2point.parameter(p_ind));
            end
            fprintf('\n');
        end
        [C1funcs,C1branches,suc]=C1branch_from_C2point(funcs,C2point,branch.parameter.free,...
            'codim2',C2type,'codim1',entry.name,'print',options.print-1,pass_on{:});
        if options.print>0
            fprintf(' emerge(s) %d x %s, ',length(C1branches),entry.name);
        end
        if isnumeric(entry.selection)
            C1branches=C1branches(entry.selection);
            if options.print>0
                for k=1:length(entry.selection)
                    fprintf(' %d ',entry.selection(k));
                end
                fprintf('selected.\n');
            end                
        end
        for k=length(C1branches):-1:1
            C1info=[C1info,struct('C1type',entry.name,'C2type',C2type,...
                'funcs',C1funcs,'branch',C1branches(k),...
                'reverse',entry.reverse,...
                'success',all(suc(:,k)))]; %#ok<AGROW>
        end
    end
end
end
