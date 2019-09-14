function [r,J]=sys_cond_extradelay(p,active,itau,iextra)
%% fix extra delays to be equal to original delays if active
% $Id: sys_cond_extradelay.m 370 2019-08-27 00:10:48Z jansieber $
%% 
[dum,isel]=ismember(active,itau); %#ok<ASGLU>
isel=isel>0;
itau=itau(isel);
iextra=iextra(isel);
ntau=length(itau);
J=repmat(p_axpy(0,p,[]),ntau,1);
if length(p)>1
    J=repmat(J,length(p),1);
end
r=p(1).parameter(itau)-p(1).parameter(iextra);
r=r(:);
for i=1:ntau
    J(i).parameter(itau(i))=1;
    J(i).parameter(iextra(i))=-1;
end
end
  
