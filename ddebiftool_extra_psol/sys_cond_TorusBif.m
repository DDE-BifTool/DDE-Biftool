function [r,J]=sys_cond_TorusBif(point,pref,ind,nonsquare)
%% additional conditions for torus and period doubling bifurcation
%
% $Id: sys_cond_TorusBif.m 369 2019-08-27 00:07:02Z jansieber $
%
Jtemplate=p_axpy(0,point,[]);
irgx=1:ind.dim;
irgu=ind.dim+(1:ind.dim);
irgv=ind.dim*2+(1:ind.dim);
irguv=[irgu,irgv];
%% keep Floquet mode at unit length
uvpoint=Jtemplate;
uvpoint.profile(irguv,:)=point.profile(irguv,:);
[utuvtv,W1]=dde_coll_profile_dot(uvpoint,uvpoint);
utuvtvJ=Jtemplate;
utuvtvJ.profile(:)=2*W1*uvpoint.profile(:);
utuvtvres=utuvtv-1;
%% ensure that real and imaginary part of Floquet mode are orthogonal
% to reference
vupoint=pref;
vupoint.profile(irgx,:)=0;
vupoint.profile([irgu,irgv],:)=vupoint.profile([irgv,irgu],:);
[utvres,Worth]=dde_coll_profile_dot(uvpoint,vupoint);
utvJ=Jtemplate;
utvJ.profile(:)=Worth*vupoint.profile(:);
%% assemble both to return array
if ~nonsquare
    r=[utuvtvres;utvres];
    J=[utuvtvJ;utvJ];
else
    r=utuvtvres;
    J=utuvtvJ;
end
end
