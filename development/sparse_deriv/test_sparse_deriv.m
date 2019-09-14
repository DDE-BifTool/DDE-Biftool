%% test Psolfold
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities');
format compact
format short g
load('minimal_demo_parnames.mat')
load('minimal_demo_stst_psol_results.mat')
%%
ipf=find(nunst_per==0,1,'first');
pf0=per_orb.point(ipf);
free_par=[];
x0=dde_x_from_point(pf0,free_par);
fcomp=@(x,str)dde_coll_jac(funcs,dde_point_from_x(x,pf0,free_par),free_par,'output',str);
f=@(x)fcomp(x,'res');
mfroms=@(s)[getfield(s,'profile'),getfield(s,'period'),...
    getfield(s,'parameter')];
Jf=@(x)mfroms(fcomp(x,'J'));
r=f(x0);
J=Jf(x0);
%%
Jp=dde_jacobian1pattern(size(J,2),J);
for i=1:max(Jp.devnumber)
    cols=Jp.i_sorted(Jp.devnumber==i,2);
    nr=find(J(:,cols));
    n2=sort(nr(:));
    j=find(diff(n2)==0,1,'first');
    fprintf('clr=%d\n',i);
    if ~isempty(j)
        fprintf('ind=%d fails\n',j);
        fail=true;
        break
    end
end
%%
J1=dde_sparsejacobian(f,x0,'pat',Jp);
Jp2=dde_jacobian2pattern(size(J,2),J);
J2=dde_sparsejacobian2(f,x0,'pattern',Jp2);
%%
free_par=[ind.b,ind.tau];
x0=dde_x_from_point(pf0,free_par);
fcomp=@(x,str)dde_coll_jac(funcs,dde_point_from_x(x,pf0,free_par),free_par,'output',str);
f=@(x)fcomp(x,'res');
mfroms=@(s)[getfield(s,'profile'),getfield(s,'period'),...
    getfield(s,'parameter')];
Jf=@(x)mfroms(fcomp(x,'J'));
nx=length(x0);
J2n=cell(1,nx);
hdev=1e-6;
for i=1:nx
    ev=zeros(nx,1);
    ev(i)=1;
    J2n{i}=(Jf(x0+ev*hdev)-Jf(x0-ev*hdev))/(2*hdev);
end
ircc=[];
J2nvals=[];
for i=1:nx
    [ir2,ic2,Jv2]=find(J2n{i});
    ircc=[ircc;ir2,ic2,i+0*ic2];
    J2nvals=[J2nvals;Jv2];
end
J1=dde_coll_dirderi(funcs,pf0,free_par,'order',1,'output','J');
norm(J1-Jf(x0),'inf')
[ircc,reo]=sortrows(ircc,[1,2,3]);
J2nvals=J2nvals(reo);
[r,J1]=dde_coll_dirderi(funcs,pf0,free_par);
J2=dde_coll_dirderi(funcs,pf0,free_par,'order',2,'output','J');
tfuncs=set_symfuncs(@sym_minimal_demo,'sys_tau',@()ind.tau,'p_vectorized',false);
J3=dde_coll_dirderi(tfuncs,pf0,free_par,'order',2,'output','J');
%%
fs=@(x)[x(1,:).*x(2,:);x(2,:).*x(3,:);x(3,:).*x(4,:)];
x1=ones(4,1);
J0=sparse(ScJacobian(fs,x1));
Jp1=dde_jacobian1pattern(4,J0);
Jp2=dde_jacobian2pattern(4,J0);
%%
dfshape=dde_coll_dirderi(funcs,pf0,[],free_par,'shape',true);
Jp=dde_jacobianpattern(dfshape);
dffunc=@(x,ydev)dde_coll_dirderi(funcs,dde_point_from_x(x,pf0,free_par),ydev,free_par);
x0=dde_x_from_point(pf0,free_par);
Jdir=@(x)dde_sparsejacobian(dffunc,x,'pattern',[],'vectorized',true,'isdfx',true);
