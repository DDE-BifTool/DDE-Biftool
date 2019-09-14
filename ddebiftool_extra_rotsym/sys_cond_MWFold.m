function [r,J]=sys_cond_MWFold(point,pref,funcs)
%% constraints used for extended DDE in fold continuation for relative periodic orbits
%
% $Id: sys_cond_MWFold.m 369 2019-08-27 00:07:02Z jansieber $
%
%% duplicate rotational phase condition (last user condition)
dim=size(point.profile,1)/2;
irgx=1:dim;
irgv=dim+(1:dim);
A=funcs.rotation;
Apref=pref;
Apref.profile(irgx,:)=0;
Apref.profile(irgv,:)=A*pref.profile(irgx,:);
[r,J]=p_dot(point,Apref,'free_par_ind',[],'period',false);
end
