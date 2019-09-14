function [col,l]=br_measr(branch,m)

% function [col,lengths]=br_measr(branch,measure)
% INPUT:
%	branch given branch
%	measure given measure
% OUTPUT:
%	col column (or matrix) containing measures along branch
%	lengths contains filled in length of each row 

% (c) DDE-BIFTOOL v. 1.00, 11/03/2000

ll=length(branch.point);

for i=1:ll
  r=p_measur(branch.point(i),m);
  l(i)=length(r);
  if size(r,1)==1
    col(i,1:l(i))=r;
  else 
    col(i,1:l(i))=r';
  end;
end;

return;
