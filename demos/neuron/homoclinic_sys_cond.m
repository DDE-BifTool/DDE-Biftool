function [res,p]=homoclinic_sys_cond(point,period)
res(1)=point.period-period;
p=p_axpy(0,point,[]);
p.period=1;
end
