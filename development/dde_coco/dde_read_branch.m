function branch=dde_read_branch(run,varargin)
%% extract DDE-Biftool branch from bifurcation diagram and one solution
% $Id: dde_read_branch.m 345 2019-05-12 00:11:15Z jansieber $
%%
default={'id',[],'label',1};
options=dde_set_options(default,varargin,'pass_on');
if isempty(options.id)
    d=coco_read_solution(options.id,run,options.label(1),'data');
    fid=d{1};
else
    fid=options.id;
end
data=coco_read_solution(fid,run,options.label(1),'data');
bd=coco_bd_read(run);
pts=coco_bd_col(bd,coco_get_id(data.id,'point'));
stid=coco_get_id(options.id,'stability');
if ismember(stid,bd(1,:))
    stab=coco_bd_col(bd,stid);
    data.info.point=cellfun(@(p,s)setfield(p,'stability',s),pts,stab)';
else
    data.info.point=[pts{:}];
end    
branch=data.info;
if isfield(branch,'tangent')
    branch=rmfield(branch,'tangent');
end
end