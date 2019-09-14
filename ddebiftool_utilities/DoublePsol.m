%% Branch off at period doubling (branch psolbranch, point number ind)
%%
function [per2branch,suc]=DoublePsol(funcs,psolbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |psolbranch|: branch of |'psol'| periodic orbits from which one wants to branch off
% * |ind|: index in |'point'| field of |psolbranch| which is close to
% period doubling where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'radius'|: initial deviation along period-doubling eigenvector
%
% All other name-value pairs are passed on to output branch.
%% Outputs
%
% * |per2branch|: branch of periodic orbits with desired settings, double period,
% and two initial corrected points
% * |suc|: flag indicating success
%
% $Id: DoublePsol.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'radius',0.01};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% create branch per2 of periodic solutions starting from an
% approximate period doubling num on a branch br of periodic orbits
psolbranch=replace_branch_pars(psolbranch,[],pass_on);
point=psolbranch.point(ind);
[per2,tangent]=dde_psol_from_psol(point,'funcs',funcs,...
    'method',psolbranch.method.stability,'radius',options.radius,pass_on{:});
free_pars=psolbranch.parameter.free;
mthp=psolbranch.method.point;
[psol1,suc]=p_correc(funcs,per2,free_pars,tangent,mthp,0,per2);
if suc==0
    per2branch=[];
    return
end
psol2ini=p_axpy(-2*options.radius,tangent,per2);
[psol2,suc]=p_correc(funcs,psol2ini,free_pars,tangent,mthp,0,psol2ini);
if suc==0
    per2branch=[];
    return
end
per2branch=setfield(psolbranch,'point',[psol2,psol1]); %#ok<SFLD>
end
