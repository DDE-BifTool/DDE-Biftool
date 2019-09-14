%% Initialize continuation branch intersecting current branch in branch point
function [branch1,suc,branch]=SetupBranchSwitch(funcs,branch,inds,varargin)
%% Input
%
% * |funcs|: functions used for DDE
% * |branch|: branch along which bifurcation was discovered
% * |inds|: indices bracketing the branch point
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of continuation parameters
%  (replacing free pars of branch)
% * |'correc'| (logical, default true): apply correction to first points on
% bifurcating branch
% * |'step'| (real, default |1e-3|): size of initial step
%
% All other named arguments are passed on to fields of |bifbranch|
%% Outputs
% 
% * |branch1|: secondary branch with 3 points (incl BP)
% * |suc|: flag whether corection was successful
% * |branch|: primary branch with BP included
%
% Parameter limits for bifbranch etc are inherited from branch, u

%
% $Id: SetupBranchSwitch.m 308 2018-10-28 15:08:12Z jansieber $
%
%% process options
default={'contpar',[],'correc',true,'step',1e-2,'hdev',1e-4};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% initialize branch of bifurcations (bifbranch)
branch=replace_branch_pars(branch,options.contpar,pass_on);
point=branch.point(inds);
branch1=[];
%% determine branch point with higher accuracy
[bp,suc,blin]=br_bp_solve(funcs,branch,inds(1),'hdev',options.hdev);
if ~suc
    return
end
%% insert point into branch
branch.point=[branch.point(1:min(inds)),bp,branch.point(max(inds):end)];
%% apply Lyapunov Schmidt reduction in bp
mth=branch.method.point;
free_pars=branch.parameter.free;
%% switch off mesh adaptation during branching
if isfield(mth,'adapt_mesh_before_correct')
    mth.adapt_mesh_before_correct=0;
end
if isfield(mth,'adapt_mesh_after_correct')
    mth.adapt_mesh_after_correct=0;
end
xfp=@(p)dde_x_from_point(p,free_pars);
pfx=@(x)dde_point_from_x(x,bp,free_pars);
%% numerical solution vector and null vectors
u0=xfp(bp);
assert(length(blin.v0)==2);
v0=xfp(blin.v0);
assert(size(blin.w0,2)==1);
w0=blin.w0;
secant=xfp(point(2))-xfp(point(1));
%% Jacobian of determining system
J=@(x)p_correc_rhs(funcs,mth,bp,free_pars,'x',x,'output','J',pass_on{:});
[t1,t2]=dde_bp_branches(J,u0,v0,w0,secant,options.hdev); %#ok<ASGLU>
branch1=rmfield(branch,'point');
s=options.step;
branch1.point=[pfx(u0-t1*s),bp,pfx(u0+t1*s)];
%% correct and add 2nd point if desired
if options.correc
    for i=[1,3]
        [p0,suc]=p_correc(funcs,branch1.point(i),free_pars,pfx(t1),...
            mth,0,branch1.point(i));
        if ~suc
            warning('SetupBranchSwitch:correc',...
                'SetupBranchSwitch: correction of point %d failed',i);
        else
            branch1.point(i)=p0;
        end
    end
end
end
