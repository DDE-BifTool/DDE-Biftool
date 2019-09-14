function p=poly_elg(m,c)
%% return value(s) of Lagrange interpolation polynomial at c in [0,1]
% function p=poly_elg(m,c);
% INPUT:
%	m lagrange degree
% 	c evaluation point in [0,1] 
% OUTPUT:
%	p values of lagrange polynomials through 0:1/m:1 at c

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000
%
% $Id: poly_elg.m 296 2018-09-24 21:36:56Z jansieber $
%
%% changed by js for vectorisation (exploited in psol_eva, psol_jac)
p=NaN(m+1,length(c));
if m==1
  p(1,:)=1-c;
  p(2,:)=c;
elseif m==2
  c1=2*(c-0.5);
  c2=c-1;
  p(1,:)=c1.*c2;
  p(2,:)=-4*c.*c2;
  p(3,:)=c.*c1;
elseif m==3
  c1=c-0.333333333333333333;
  c2=c-0.666666666666666667;
  c3=c-1;
  c12=4.5*c1.*c2;
  c03=13.5*c.*c3;
  p(1,:)=-c12.*c3;
  p(2,:)=c03.*c2;
  p(3,:)=-c1.*c03;
  p(4,:)=c.*c12;
elseif m==4
  p(1,:)=10.666666666666666667*(c-0.25).*(c-0.5).*(c-0.75).*(c-1);
  p(2,:)=-42.666666666666666667*c.*(c-0.5).*(c-0.75).*(c-1);
  p(3,:)=64*c.*(c-0.25).*(c-0.75).*(c-1);
  p(4,:)=-42.666666666666666667*c.*(c-0.25).*(c-0.5).*(c-1);
  p(5,:)=10.666666666666666667*c.*(c-0.25).*(c-0.5).*(c-0.75);
else
  t=0:1/m:1;
  for j=1:m+1
    p(j,:)=ones(1,length(c));
    for k=1:m+1
      if k~=j
        p(j,:)=p(j,:).*(c-t(k))/(t(j)-t(k));
      end;
    end;
  end;
end
%% (js) for compatibility with original code
if length(c)==1
    p=p';
end
return;
