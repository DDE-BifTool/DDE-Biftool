function [funcs,method,free_par]=dde_point_delay_zero_prep(funcs,method,free_par,d_nr,tz)
%% adjust Newton iteration parameters for finding solution with zero delay 
% for collocation mesh (state-dependent delays)
%
% Inputs:
%
% * funcs: rhs and delays
% * method: Newton iteration parameters
% * free_par: free parameter indices along branch
% * d_nr: number of delay becoming negative
% * t_z: time along mesh where sign change of delay should happen
%
% $Id: dde_point_delay_zero_prep.m 326 2019-02-03 14:54:42Z mmbosschaert $
%%
method.extra_condition=1;
if isfield(method,'phase_condition')
    method.phase_condition=0;
end
type=method.delay_zero_prep;
funcs.sys_cond=@(pt)extcond(pt,type,funcs.sys_cond,funcs,free_par,d_nr,tz);
end
function [r,J]=extcond(pt,type,usercond,funcs,free_par,d_nr,tz)
[r1,J1]=usercond(pt);
[r2,J2]=dde_apply({'dde_',type,'_delay_zero_cond'},funcs,pt,free_par,d_nr,tz);
r=[r1(:);r2(:)];
if isempty(J1)
  J=[J2(:)];
else
  J=[J1(:);J2(:)];
end
end
