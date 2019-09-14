function [resi,condi]=hom_cond(point)

resi=point.parameter(7)-point.period;

condi=p_axpy(0,point,[]);
condi.period=-1;
condi.parameter(7)=1;
return;
