function branch=br_stabl(funcs,branch,skip,recompute)
%% compute stability information along branch
% function st_branch=br_stabl(funcs,branch,skip,recompute)
% INPUT:
%   funcs problem function
%	branch 
%	skip number of points to skip between stability computations
%	recompute if nonzero recompute stability info already present
% OUTPUT:
%	st_branch branch with stability information

% (c) DDE-BIFTOOL v. 1.02, 02/11/2000
%
% $Id: br_stabl.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
ll=length(branch.point);

if ll<1 
  err=ll;
  error('BR_STABL: branch is empty: %d points!',err);
end;

if ~isfield(branch.point(1),'stability')
  branch.point(1).stability=[];
end;

for i=1:skip+1:ll-1
  if isempty(branch.point(i).stability) || recompute
    branch.point(i).stability=p_stabil(funcs,branch.point(i),branch.method.stability);
  end;
end;

if isempty(branch.point(ll).stability) || recompute
  branch.point(ll).stability=p_stabil(funcs,branch.point(ll),branch.method.stability);
end;

return;
