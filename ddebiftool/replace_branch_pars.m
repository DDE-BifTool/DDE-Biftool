function branch=replace_branch_pars(branch,contpar,pass_on)
%% pass on optional arguments to branch structure
%
% branch=replace_branch_pars(branch,contpar,pass_on)
% 
% input
% branch: given branch structure to be amended
% contpar: continuation parameter, prepended to free parameters already
%          present in branch (if length==1) or replacing
%          branch.parameter.free
% pass_on: cell array containing name-value pairs for fields of
%          branch.method.continuation, branch.method.point, branch.method.stability and
%          branch.parameter to be replaced
%
% output: branch with amended fields
%
% $Id: replace_branch_pars.m 339 2019-05-09 19:47:01Z jansieber $
%
%%
branch.method=replace_method_pars(branch.method,pass_on{:});
if isempty(contpar)
    % if no continuation parameters are given use free parameters of branch
    branch.parameter.free=branch.parameter.free;
else
    branch.parameter.free=contpar(:)';
end
branch.parameter=dde_set_options(branch.parameter,pass_on,'pass_on');
end
