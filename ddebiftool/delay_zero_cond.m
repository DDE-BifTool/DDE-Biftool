function [r,Jp]=delay_zero_cond(funcs,point,d_nr,free_par_ind,tz)
%% residual & derivative wrt x, T and parameter for tau_d_nr(tz)=tau'_dnr(tz)=0
% formulated in sys_cond format
%
% inputs
% funcs: system functions
% point: point (type psol,stst,hopf,fold) for which delay=0 is to be
% determined
% tz: time (for psol) at which delay is supposed to be zero (for psol)
% d_nr: delay number
% free_par_ind: indices of free parameters
%
% outputs
% r residual: tau for stst,hopf,fold (1d), tau(tz) and tau'(tz) for psol
% J Jacobian (struct of type point): single row for stst,hopf,fold, two
% rows for psol
%
% $Id: delay_zero_cond.m 296 2018-09-24 21:36:56Z jansieber $
%
ntau=funcs.sys_ntau();
if strcmp(point.kind,'stst') || strcmp(point.kind,'hopf') || strcmp(point.kind,'fold')
    xx=point.x(:,ones(ntau+1,1));
    r=funcs.sys_tau(d_nr,xx,point.parameter);
    Jp=p_axpy(0,point,[]);
    for j=0:ntau
        Jp.x=Jp.x+funcs.sys_dtau(d_nr,xx,point.parameter,j,[])';
    end
    for k=1:length(free_par_ind)
        Jp.parameter(free_par_ind(k))=funcs.sys_dtau(d_nr,xx,point.parameter,[],free_par_ind(k));
    end
elseif strcmp(point.kind,'psol')
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
    [J,r]=psol_jac(tfuncs,t,point.period,point.profile,point.mesh,point.degree,...
        point.parameter,free_par_ind,[],...
        'c_is_tvals',true,'bc',false,'Dtmat',zeros(1,size(point.profile,1)));
    r=[r(1);(r(3)-r(2))/2/h];
    J=[J(1,:);(J(3,:)-J(2,:))/2/h];
    Jp=p_axpy(0,point,[]);
    Jp=[Jp;Jp];
    nx=numel(point.profile);
    np=length(free_par_ind);
    for k=1:2
        Jp(k).profile(:)=J(k,1:nx);
        Jp(k).period=J(k,nx+1);
        Jp(k).parameter(free_par_ind)=J(k,nx+1+(1:np));
    end
else
    error('delay_zero_cond:kind','delay_zer_cond: kind %s not implemented',...
        point.kind);
end
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
