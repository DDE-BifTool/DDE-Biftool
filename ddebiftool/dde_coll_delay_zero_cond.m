function [r,Jp]=dde_coll_delay_zero_cond(funcs,point,free_par_ind,d_nr,tz)
%% residual & derivative wrt x, T and parameter for tau_d_nr(tz)=tau'_dnr(tz)=0
% formulated in sys_cond format
%
% inputs
% funcs: system functions
% point: point (type psol) for which delay=0 is to be
% determined
% tz: time for psol at which delay is supposed to be zero (for psol)
% d_nr: delay number
% free_par_ind: indices of free parameters
%
% outputs
% r residual: tau for stst,hopf,fold (1d), tau(tz) and tau'(tz) for psol
% J Jacobian (struct of type point): single row for stst,hopf,fold, two
% rows for psol
%
% $Id: dde_coll_delay_zero_cond.m 366 2019-07-14 21:52:54Z jansieber $
%
tfuncs=funcs;
tfuncs.sys_ntau=@()d_nr;
tfuncs.sys_rhs=@(xx,p)funcs.sys_tau(d_nr,xx(:,1:d_nr+1,:),p);
f=@(xx,p,nx,np,v)funcs.sys_dtau(d_nr,xx(:,1:d_nr+1,:),p,nx,np);
tfuncs.sys_deri=@(xx,p,nx,np,v)sys_deri(xx,p,nx,np,v,f);
it0=find(point.mesh<=tz,1,'last');
if point.mesh(it0)==tz
    itlow=it0-1;
else
    itlow=it0;
end
itup=itlow+1;
if itlow<1
    tlow=point.mesh(end-1)-1;
else
    tlow=point.mesh(itlow);
end
if itup>length(point.mesh)
    tup=point.mesh(2)+1;
else
    tup=point.mesh(itup);
end
h=min([tz-tlow,tup-tz]);
t=[tz,tz-h,tz+h];
[J,r]=dde_coll_jac(tfuncs,point,free_par_ind,'c',t,...
    'c_is_tvals',true,'Dtmat',zeros(1,size(point.profile,1)),'wrapJ',true);
dh=@(M)(M(3,:)-M(2,:))/2/h;
rs=@(M)reshape(M,size(point.profile));
r=[r(1);dh(r)];
Jp=repmat(p_axpy(0,point,[]),2,1);
Jp(1).profile=rs(J.profile(1,:));
Jp(1).period=J.period(1);
Jp(1).parameter(free_par_ind)=J.parameter(1,:);
Jp(2).profile=rs(dh(J.profile));
Jp(2).period=dh(J.period);
Jp(2).parameter(free_par_ind)=dh(J.parameter);
end
%%
function J=sys_deri(xx,p,nx,np,v,f)
J=f(xx,p,nx,np,v);
if ~isempty(nx) && isempty(np)
    J=reshape(J,[1,size(xx,1),size(xx,3)]);
elseif isempty(nx) && ~isempty(np)
    J=reshape(J,[1,1,size(xx,3)]);
else
    error('delay_zero_cond:order',...
        'delay_zero_cond: higher-order derivatives not implemented')
end
end
