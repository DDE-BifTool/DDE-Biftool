function [res,J]=rot_cond(point,funcs,pref)
%% extra phase condition needed to fix rotation phase
%
% $Id: rot_cond.m 357 2019-06-30 23:59:31Z jansieber $
%
A=funcs.rotation;
usercond=funcs.orig_cond;
userpoint=point;
userpoint.parameter=userpoint.parameter(1:end-1);
if funcs.orig_cond_reference
    userref=pref;
    userref.parameter=userref.parameter(1:end-1);
    [userres,userJ]=usercond(userpoint,userref);
else
    [userres,userJ]=usercond(userpoint);
end
if isfield(point,'x')
    resphas=pref.x'*A*point.x;
    Jphas=p_axpy(0,point,[]);
    Jphas.x=A'*pref.x;
elseif isfield(point,'profile')
    pref.profile=A*pref.profile;
    pref.period=0;
    [resphas,Jphas]=p_dot(point,pref,'free_par_ind',[]);
else
    error('rot_cond:type','rot_cond: type %s not supported',point.kind);
end
res=[userres(:);resphas];
if ~isempty(userJ) % fix for octave
    J=[userJ(:);Jphas];
else
    J=Jphas;
end
end
