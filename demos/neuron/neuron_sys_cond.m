function [res,p]=neuron_sys_cond(point)
res(1)=point.x(1);
p=p_axpy(0,point,[]);
p.x(1)=1;

return;
