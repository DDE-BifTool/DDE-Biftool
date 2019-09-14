function br_bifplot(branch,x_m,y_m,method)
%% plot branch
% function br_bifplot(branch,x_measure,y_measure,line_type)
% INPUT:
%	branch: branch of points
%	x_measure: scalar measure for x-coordinate
%	y_measure: scalar measure for y_coordinate
%	varargin: list of type 'bifurcation', 'line type', 'bifurcation',
%	'line type'
%
%  $Id: br_bifplot.m 309 2018-10-28 19:02:42Z jansieber $
%%

FPI = br_getflags(branch);
[totalbifnum,~] = size(FPI);
if isempty(FPI)
   fprintf('BR_BIFPLOT: no bifurcations present in branch!\n');
   return
end

%% Check all bifurcation types
for i = 1:totalbifnum
   bifcount = length(FPI(i,:)); 
   if bifcount <= 0
      continue
   end
   linetype = method.(num2bif(i));
   if isempty(linetype)
       continue
   end
   % Plot all points
   for j = 1:bifcount
      pointno = FPI(i,j);
      if pointno <= 0
          break
      end
      point = branch.point(pointno);
      y = p_measur(point,y_m);
      x = p_measur(point,x_m);
      plot(x,y,linetype);
   end
end
end

