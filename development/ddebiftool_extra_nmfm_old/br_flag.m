function flaggedbranch = br_flag(branch)
%% set empty flag on all points of branch
% function flaggedbranch = br_flag(branch)
% Purpose:
%   Sets flag = '' on all points of the branch.
% INPUT:
%	branch 
% OUTPUT:
%	flaggedbranch
%
% $Id: br_flag.m 65 2014-12-23 22:37:36Z jan.sieber $
%

flaggedbranch = branch;
ll=length(branch.point);

if ll<1
   error('BR_FLAG: branch is empty!');
end

for i=1:ll
   if ~isfield(flaggedbranch.point(i),'flag') || isempty(flaggedbranch.point(i).flag)
      flaggedbranch.point(i).flag = '';
   end
   if ~isfield(flaggedbranch.point(i),'nmfm')
      flaggedbranch.point(i).nmfm = [];
   end
   if ~isfield(flaggedbranch.point(i),'nvec')
      flaggedbranch.point(i).nvec = [];
   end
end
end
