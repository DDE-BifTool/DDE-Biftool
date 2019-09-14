function [px]=hcli_eva(profile1,t,x,m);

% function [px]=hcli_eva(profile1,t,x,m);
% INPUT:
%	profile1 profile on mesh t
%	t representation points in [0,1]^(m*l+1)
%	x point(s) where to evaluate
%	m order polynomials
% OUTPUT: 
%	px value of profile at x

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

for i=1:length(x)
  if x(i)>t(end)
     px(:,i)=profile1(:,end);
  elseif x(i)<t(1)
    px(:,i)=profile1(:,1);
  else
    c=x(i);
    index=length(t)-m;
    while t(index)>c
      index=index-m;
    end;
    Pa=poly_elg(m,(c-t(index))/(t(index+m)-t(index)));
    px(:,i)=profile1(:,index:index+m)*Pa';
  end;
end;

return;
