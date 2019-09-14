function initbranch = br_bifinit(branch)
%% Add flag, nmfm and nvec on all points of Hopf branch.
% function initbranch = br_bifinit(branch)
% Purpose:
%   Adds flag, nmfm and nvec (hopf only) on all points of the branch.
% INPUT:
%	branch 
% OUTPUT:
%	initbranch
%
% $Id: br_bifinit.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
initbranch = branch;
ll=length(branch.point);
kind = branch.point(1).kind;

if ll<1
   error('BR_FLAG: branch is empty!');
end

for i=1:ll
   if ~isfield(initbranch.point(i),'flag') || isempty(initbranch.point(i).flag)
      initbranch.point(i).flag = '';
   end
   if ~isfield(initbranch.point(i),'nmfm')
      initbranch.point(i).nmfm = struct();
   end
   if strcmp(kind,'hopf') && ~isfield(initbranch.point(i),'nvec')
      initbranch.point(i).nvec = struct();
   end
end
end
