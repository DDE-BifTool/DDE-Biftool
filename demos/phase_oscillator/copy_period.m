%% Use period as parameter (which can be fixed or released)
%
% $Id: copy_period.m 126 2016-09-05 22:37:08Z jansieber $
%
%%
function [res,J]=copy_period(pt,indpar)
res=pt.period-pt.parameter(indpar);
J=p_axpy(0,pt,[]);
J.period=1;
J.parameter(indpar)=-1;
end