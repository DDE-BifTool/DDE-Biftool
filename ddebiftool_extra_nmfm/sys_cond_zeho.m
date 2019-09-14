function [res,J]=sys_cond_zeho(point,get_comp,usercond)
%% extra conditions added for nonlinear system determining zeho point
%
% 2 conditions are needed: $q_0^T q_0 = 1$, $q_1^T q_1 = 1$ where $q_0$ is
% the nullvector of $D(0)$ (the linearized system) and $q_1$ is the
% nullvector of $D_om$ (see sys_rhs_zeho.m)
%
% $Id: sys_cond_zeho.m 327 2019-02-05 12:18:54Z mmbosschaert $
%%
if ~strcmp(point.kind,'hopf') && ...
        ~strcmp(point.kind,'fold') && ...
        ~strcmp(point.kind,'stst')
    error('sys_cond_zeho:kind',['sys_cond_zeho was called for %s, but ',...
        'can only be called for stst,fold or hopf'],point.kind);
end
%% user conditions
userpoint=get_comp(point,'solution');
[userres,userJ]=usercond(userpoint);
%% append artificial components to user conditions as zeros
ind_x=get_comp(point,'ind_x');
ind_q0=get_comp(point,'ind_q0');
ind_rq1=get_comp(point,'ind_rq1');
ind_iq1=get_comp(point,'ind_iq1');
ind_p=get_comp(point,'ind_parameter');
uzhJ=repmat(point,length(userJ),1);
for i=length(userJ):-1:1
    uzhJ(i)=p_paxpy(0,point,[]);
    uzhJ(i).x(ind_x)=userJ(i).x;
    uzhJ(i).parameter(ind_p)=userJ(i).parameter;
end
%% zeho extra conditions
q0=get_comp(point,'q0');
q1=get_comp(point,'q1');
zhres=[q0'*q0-1;q1'*q1-1;0];
zhJ=repmat(p_axpy(0,point,[]),3,1);
% derivaties for <q0,q0> -1 and <q1,q1>
zhJ(1).x(ind_q0)=2*q0;
zhJ(2).x(ind_rq1)=2*real(q1);
zhJ(2).x(ind_iq1)=2*imag(q1);
zhJ(3).x(ind_rq1)=-imag(q1);
zhJ(3).x(ind_iq1)=real(q1);

%% assemble conditions
res=[userres;zhres];
J=[uzhJ;zhJ];
end
