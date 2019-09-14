function [r,J]=sys_cond_POfold(point,ind)
%% constraints used for extended DDE in periodic fold continuation
%
% $Id: sys_cond_POfold.m 369 2019-08-27 00:07:02Z jansieber $
%
%% obtain condition that nullvector has length 1
vpoint=p_axpy(0,point,[]);
vpoint.profile(ind.dim+1:end,:)=point.profile(ind.dim+1:end,:);
vpoint.parameter=point.parameter;
[vtv,vtvJ]=p_dot(vpoint,vpoint,'free_par_ind',[ind.beta,ind.nullparind(:,2)']);
vtvJ=p_axpy(2,vtvJ,[]);
vtvres=vtv-1;
%% obtain condition <x',v>=0
p1=setfield(point,'profile',point.profile(ind.dim+(1:ind.dim),:));  %#ok<*SFLD>
p2=setfield(point,'profile',point.profile(1:ind.dim,:));
[xdtvres,pW]=dde_coll_profile_dot(p1,p2,'derivatives',[0,1]);
xdtvJ=p_axpy(0,point,[]);
xdtvJ.profile=cat(1,reshape(pW'*p1.profile(:),ind.dim,[]),...
    reshape(pW*p2.profile(:),ind.dim,[]));
%% assemble residuals and Jacobians
r=[vtvres;xdtvres];
J=[vtvJ;xdtvJ];
end
