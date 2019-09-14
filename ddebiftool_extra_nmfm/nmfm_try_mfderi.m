function [success,y]=nmfm_try_mfderi(funcs,point,devs)

ind_tau=funcs.sys_tau();

xx=point.x;
xx=xx(:,ones(1,length(ind_tau)+1));
tau0=[0,point.parameter(ind_tau)];
success=true;
devvals={};
for i=length(devs):-1:1
    devvals{i}=nmfm_dev_call(devs(i),-tau0);
    devvals{i}=devvals{i}(:);
end
    y=funcs.sys_mfderi(xx,point.parameter,devvals{:});
end