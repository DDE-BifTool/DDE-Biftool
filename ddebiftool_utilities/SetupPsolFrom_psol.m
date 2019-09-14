%% Start branch of periodic orbits from bifurcation of periodic orbit branch
%%
function [per,suc]=SetupPsolFrom_psol(bfuncs,branch,ind,varargin)
%% Inputs
% 
% * |bfuncs|: structure with extended functions or functions provided by
% user
% * |branch|: branch of |'psol'| solutions that are bifurcations of psols
% * |ind|: index in |'point'| field of |branch| where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% (general)
%
% * |'contpar'|: index of continuation parameters if different from
% |'branch'|, (this should be set manually when the original branch was a
% two-parameter branch)
% * |'corpar'|: parameters left free for initial correction (if different
% from |contpars|)
%
% All other optional inputs are passed on to fields of out branch |per|
%% Outputs
%
% * |per|: branch of periodic orbits with desired settings and two initial
% corrected points
% * |suc|: flag indicating success
% 
% $Id: SetupPsolFrom_psol.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'contpar',[],'corpar',[],'step',1e-3,'correction',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if ~isfield(bfuncs,'userfuncs')
    %% this is a branch of periodic orbits
    [per,suc]=ChangeBranchParameters(bfuncs,branch,ind,'contpar',options.contpar,...
        'corpar',options.corpar,pass_on{:});
else
    funcs=bfuncs.userfuncs;
    psol=bfuncs.get_comp(branch.point(ind),'solution');
    pbranch=setfield(branch,'point',psol); %#ok<SFLD>
    pbranch.method.point.preprocess='';
    pbranch.method.point.postprocess='';
    pbranch=replace_branch_pars(pbranch,options.contpar,pass_on);
    [per,suc]=correct_ini(funcs,pbranch,psol,options.contpar,options.step,...
        options.correction);
end
end
