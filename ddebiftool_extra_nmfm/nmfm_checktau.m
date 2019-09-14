%% check if parameter deviations are delays
function funcs=nmfm_checktau(funcs,point,devs)
if ~funcs.tp_del
    itau=funcs.sys_tau();
    ndim=size(point.x,1);
    for i=length(devs):-1:1
        dev0(:,i)=nmfm_dev_call(devs(i),0);
    end
    pdevs=dev0(ndim+1:end,:)';
    if any(any(pdevs(:,itau)))
        %% we need delay functions
        funcs.sys_tau=@(it,xx,p)nmfm_taufunc(0,it,xx,p,[],itau);
        maxorder=length(devs);
        for order=maxorder:-1:1
            funcs.sys_dirdtau{order}=@(it,xx,p,dxx,dp)nmfm_taufunc(order,it,xx,p,dp,itau);
        end
        funcs.sys_ntau=@()length(itau);
        funcs.tp_del=true;
        funcs.sys_dirdtau_provided=maxorder;
    end
end
end
%%