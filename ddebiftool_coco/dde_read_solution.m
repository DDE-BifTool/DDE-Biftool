function data=dde_read_solution(run,varargin)
%% wrapper for coco_read_solution, returns DDE-Biftool type branch structure
% $Id: dde_read_solution.m 346 2019-05-13 05:41:50Z jansieber $
%
%% Inputs
%
% * |run|: name of run
% * |'id'|: prefix for instance (can be empty)
% * |'label'|: if number, extract this label number, if character, extract
% all labels with this flag (e.g., 'HB' extracts all Hopf points).
%
%% Output
% data saved for this solution. chart data is put into data.info.point and
% data.info.tangent. This can be passed on to |dde_construct|.
% 
%%
default={'id',[],'label',1};
options=dde_set_options(default,varargin,'pass_on');
if ischar(options.label)
    bd=coco_bd_read(run);
    labels=coco_bd_labs(bd,options.label);
else
    labels=options.label;
end
if isempty(labels)
    data=[];
    return
end
if isempty(options.id)
    d=coco_read_solution(options.id,run,labels(1),'data');
    options.id=d{1};
end    
for i=length(labels):-1:1
    [data(i),chart]=coco_read_solution(options.id,run,labels(i),'data','chart');
    data(i).info.point=dde_point_from_x(chart.x,data(i).info.point,data(i).ipar);
    data(i).info.tangent=dde_point_from_x(chart.t,data(i).info.point,data(i).ipar);
end
end
