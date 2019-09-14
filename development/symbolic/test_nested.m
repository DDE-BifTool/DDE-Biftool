%% test symbolic toolbox for state-dependent delays
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');
%%
ntau=2;
nx=2;
maxorder=5;
xx=sym('xx',[nx,ntau+1]);
par=sym('p',[1,2]);
f=[par(1)-xx(1,ntau+1);-xx(2,ntau+1)];
for i=1:ntau
    tau(i)=xx(1,i);
end
%%
%[fstr,df,xxdev,pdev]=dde_symdericode(f,xx,par,'filename','sym_nested_rhs');
%%
%[taustr,dtau,xxtdev,ptdev]=dde_symdericode(tau,xx,par,'filename','sym_nested_tau','scalar',true);
%%
%[df,xxdev]=dde_sd_symdiff(f,tau,xx);
%%
[funcstr,derivs]=dde_sym2funcs(f,xx,par,'sd_delay',tau,'filename','sym_nested');
%% replace xx by function


t=sym('t');
h=sym('h');
clear xf
for i=nx:-1:1
    xfname{i}=['x',num2str(i)];
    vfname{i,1}=['v',num2str(i)];
    xf{i,1}=symfun([xfname{i},'(t)'],t);
    vf{i,1}=symfun([vfname{i},'(t)'],t);
    for k=1:maxorder
        vfname{i,k+1}=[vfname{i,1},'_',num2str(k)];
        vf{i,k+1}=symfun([vfname{i,k+1},'(t)'],t);
    end
end
tauf=[0;tau(:)];
ff=f;
for k=1:ntau+1
    for i=1:nx
        for l=1:k-1
            tauf(k)=subs(tauf(k),xx(i,l),xf{i}(t-tauf(l)));
        end
        ff=subs(ff,xx(i,k),xf{i}(t-tauf(k)));
    end
end
%%
fd=ff;
for i=1:nx
    fd=subs(fd,xf{i},xf{i}+h*vf{i});
end
%% remove derivatives of x (in equlibrium)
fdiff=diff(fd,h,h,h);
fdiff=subs(fdiff,h,0);
xeq=sym('xeq',[nx,1]);
fdeq=fdiff;
for i=1:nx
    fdeq=subs(fdeq,xf{i},xeq(i));
end
%% replace derivatives of v with dvk
fdeqs=fdeq;
for i=1:maxorder
    subseqs=cell(1,nx);
    for k=1:nx
        subseqs{k}=['D(',vfname{k,i},')=',vfname{k,i+1}];
    end
    fdeqs=feval(symengine,'subs',char(fdeqs),subseqs{:});
end
w_orig=warning;
warning('off','symbolic:generate:FunctionNotVerifiedToBeValid');
matlabFunction(fdeqs,'Vars',{xeq,vfname{:}},'File','nested_fcn');
warning(w_orig);
