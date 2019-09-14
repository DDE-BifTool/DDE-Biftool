function branch=br_rvers(branch)

% function t_branch=br_rvers(branch)
% INPUT:
%       branch 
% OUTPUT:
%       t_branch branch with points in reversed order

% (c) DDE-BIFTOOL v. 1.00, 11/03/2000

branch.point=branch.point(length(branch.point):-1:1);

return;

