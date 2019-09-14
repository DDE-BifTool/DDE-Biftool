function branchlist=nmfm_branch_list(type)
%% Return list of branches emerging from bifurcation
% (typically codim 1 branches emerging from codim 2 points)
%% Input
% * |type|: flag used for codim 2 point. Supported: hoho, genh, zeho. All
% others return no branches
%% Output
% * |branchlist|: array of structures with fields |name| (used to construct
% function call), |selection| (which branches to select, for example, in a
% Hopf-Hopf interaction, one wants only the first Hopf bifurcation, as the
% second one is the current branch), |reverse| (logical, a hint whether
% this bifurcation should be continued in both directions)
%
% $Id: nmfm_branch_list.m 309 2018-10-28 19:02:42Z jansieber $
%%
switch type
    case 'hoho'
        branchlist=struct('name',{'hopf','TorusBifurcation'},...
            'selection',{1,1:2},'reverse',{true,false});
    case 'genh'
        branchlist=struct('name','POfold','selection',1,'reverse',false);
    case 'zeho'
        branchlist=struct('name',{'hopf','fold','TorusBifurcation'},...
            'selection',{1,1,1},'reverse',{false,true,false});
    otherwise
        branchlist=struct('name',{},'selection',{},'reverse',{});
end
end

        
        