function [res,J]=sys_cond_BT(point,get_comp,usercond)
if ~strcmp(point.kind,'hopf') && ...
        ~strcmp(point.kind,'fold') && ...
        ~strcmp(point.kind,'stst')
    error('sys_cond_BT:kind',['sys_cond_BT was called for %s, but ',...
        'can only be called for stst,fold or hopf'],point.kind);
end
n=size(point.x,1)/3;
%% user conditions
userpoint=get_comp(point,'solution');
[userres,userJ]=usercond(userpoint);
%% append artificial components to user conditions as zeros
for i=1:length(userJ)
    userJ(i).x=[userJ.x;zeros(2*n,1)];
end
%% BT extra conditions
q0=point.x(n+(1:n));
q1=point.x(2*n+(1:n));
BTres=[q0'*q0-1;q0'*q1];
BTJ(2)=p_axpy(0,point,[]);
BTJ(1)=BTJ(2);
% derivaties for <q0,q0> -1 and <q0,q1>
BTJ(1).x(n+(1:n))=2*q0;
BTJ(2).x(n+(1:n))=q1;
BTJ(2).x(2*n+(1:n))=q0;
%% assemble conditions
res=[userres;BTres];
J=[userJ;BTJ];
end
