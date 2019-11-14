%% create mesh from coarse mesh
% if argument grid is integer then "grid" mesh points are evenly
% distributed in between (otherwise, these should be the internal mesh
% points, eg, legpts rescaled to [0,1]). 'doublecount' (false) duplicates
% boundaries, 'acc' (1) adds coarse mesh points
%
% $Id: dde_coll_meshfill.m 369 2019-08-27 00:07:02Z jansieber $
%%
function newmesh=dde_coll_meshfill(tcoarse,degree,varargin)
%% fill coarse grid with uniform values from grid
default={'acc',1,'doublecount',false,'scal',@(x)x,'grid','linear',...
    'purpose','storage'};
options=dde_set_options(default,varargin,'pass_on');
grid=dde_coll_set_grid(options.purpose,degree,'type',options.grid);
if ~options.doublecount && grid(1)==0 && grid(end)==1
    append=tcoarse(end);
    grid=grid(1:end-1);
else
    append=[];
end

tcoarse=tcoarse(:)';
grid=grid(:);
dt=diff(tcoarse);
ndt=length(dt);
og=ones(length(grid),1);
ot=ones(ndt,1);
scal=options.scal(dt);
addmesh=grid(:,ot).*scal(og,:);
repmesh=tcoarse(og,1:end-1)*options.acc+addmesh;
newmesh=[repmesh(:)',append];
end
