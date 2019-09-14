function [r,J]=sys_cond_psolfold(point,userfuncs,conversion,method)
psolfold=dde_coll_convert(point,conversion);
pfzero=p_axpy(0,psolfold,[]);
free_par=psolfold.free_par;
per=dde_psolfold_extract(psolfold,'solution');
per_var=dde_psolfold_extract(psolfold,'var');
%% incorporate user sys_cond's and their derivatives
xper=dde_x_from_point(per,free_par);
nx=length(xper);
if isfield(userfuncs,'sys_cond_deri_provided') &&userfuncs.sys_cond_deri_provided
    [r_u,J_u]=userfuncs.sys_cond([per,per_var]);
else
    [r1,J1]=userfuncs.sys_cond(per);
    n1=length(r1);
    J1x=reshape(dde_x_from_point(J1,free_par),nx,n1);
    per_dev1=p_axpy(method.hdev,per_var,per);
    [rv1,Jv1]=userfuncs.sys_cond(per_dev1); %#ok<ASGLU>
    per_dev2=p_axpy(-method.hdev,per_var,per);
    [rv2,Jv2]=userfuncs.sys_cond(per_dev2); %#ok<ASGLU>
    r_u=[r1(:),J1x'*dde_x_from_point(per_var,free_par)];
    J_e=J1(:);
    for i=1:length(J_e)
        J_e(i)=p_axpy(-1,Jv2(i),Jv1(i));
        J_e(i)=p_axpy(1/(2*method.hdev),J_e(i),[]);
    end
    J_u=[J1(:),J_e];
end
n_u=size(r_u,1);
r_user=r_u(:);
J_user=repmat(pfzero,length(r_u),1);
for i=1:length(r1)
    J_user(i).profile=  J_u(i,1).profile;
    J_user(i).parameter=J_u(i,1).parameter;
    J_user(i).period=   J_u(i,1).period;
    J_user(i+n_u).profile=       J_u(i,2).profile;
    J_user(i+n_u).parameter=     J_u(i,2).parameter;
    J_user(i+n_u).period=        J_u(i,2).period*point.period+...
                                         J_u(i,1)*per_var.period;
    J_user(i+n_u).profile_var=   J_u(i,1).profile;    
    J_user(i+n_u).parameter_var=[J_u(i,1).period/point.period,...
                                   J_u(i,1).parameter(free_par)];
end
%% extra conditions needed for psolfold
%% obtain condition <v,v>+beta^2-1=0
vpoint=per;
vpoint.profile=psolfold.profile_var;
[vtv,W]=dde_coll_profile_dot(vpoint,vpoint);
J_vtv=pfzero;
J_vtv.profile_var(:)=2*vpoint.profile(:)'*W;
J_vtv.parameter_var=2*psolfold.parameter_var;
r_vtv=vtv+sum(psolfold.parameter_var.^2)-1;
%% obtain condition <x',v>=0
[r_xdtv,Wd]=dde_coll_profile_dot(per,vpoint,'derivatives',[1,0]);
J_xdtv=pfzero;
J_xdtv.profile_var(:)=Wd'*per.profile(:);
J_xdtv.profile(:)=Wd*per_var.profile(:);
%% constrain extra delays
extra_delay=psolfold.parameter_delay;
ndel=length(extra_delay);
itau=[0,userfuncs.sys_tau()];
Jfixd=repmat(pfzero,ndel,1);
rfixd=zeros(ndel,1);
if ~userfuncs.tp_del
    rfixd(1)=-psolfold.parameter_delay(1);
    Jfixd(1).parameter_delay(1)=-1;
    for i=2:ndel
        rfixd(i)=psolfold.parameter(itau(i))-extra_delay(i);
        Jfixd(i).parameter(itau(i))=1;
        Jfixd(i).parameter_delay(i)=-1;
    end
end
r=[r_user;r_vtv;r_xdtv;rfixd];
J=[J_user(:);J_vtv;J_xdtv;Jfixd];
J=arrayfun(@dde_coll_convert,J);
end