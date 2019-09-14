function data = dde_point_update(~, data, cseg, varargin)
%% dde_point_update: set reference point
% $Id: dde_point_update.m 345 2019-05-12 00:11:15Z jansieber $
%%
% update information about current solution 
u             = cseg.src_chart.x; % Current chart
data.info.point= dde_point_from_x(u,data.info.point,data.ipar);
end
