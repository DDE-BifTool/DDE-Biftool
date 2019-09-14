function [J,res,tT,extmesh,Jstruc]=psol_jac_complete(funcs,psol,free_par,varargin)
%% optional
default={'bc',true,'ph',true,'rotationcheck',true,'step_cnd',[],'extra_cond',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[J,res,tT,extmesh,Jstruc]=psol_jac_reduced(funcs,psol,free_par,pass_on{:});
%% boundary condition:
if options.bc 
   if islogical(options.bc)
       Jbc_dx=zeros(Jstruc.cols.nx,Jstruc.cols.nx,length(Jstruc.cols.x));
       Jbc_dp=zeros(Jstruc.cols.nx,Jstruc.cols.np);
       Jbc_dT=zeros(Jstruc.cols.nx,Jstruc.cols.nT);
       Jbc_dx(:,:,1)=eye(Jstruc.cols.nx);
       Jbc_dx(:,:,end)=-eye(Jstruc.cols.nx);
       Jbc=[reshape(Jbc_dx,Jstruc.cols.nx,[]),Jbc_dT,Jbc_dp];
       nbc=Jstruc.cols.nx;
       resbc=psol.profile(:,1)-psol.profile(:,length(Jstruc.cols.x));
       if options.rotationcheck
           resbc=mod(resbc+pi,2*pi)-pi;
       end
   elseif isa(options.bc,'function_handle')
        [resbc,Jbc]=options.bc(psol.profile(:,[1,end]),psol.parameter,psol.period);
        nbc=size(resbc,1);
        Jbc_dx(:,:,1)=Jbc.x1;
        Jbc_dx(:,:,end)=Jbc.x2;
        Jbc=reshape(Jbc_dx,nbc,[]);
        if Jstruc.cols.nT>0
            Jbc=[Jbc,Jbc.T];
        end
        if Jstruc.cols.np>0
            Jbc=[Jbc,Jbc.p];
        end
   end
   res=[res;resbc(:)];
   Jstruc.rows.bc=size(J,1)+(1:nbc);
   J=[J;Jbc];
else
    nbc=0;
    Jstruc.rows.bc=size(J,1)+(1:nbc);
end
%% phase condition:
if options.ph 
    p0=p_axpy(0,psol,[]);
    [resph,p_Jph_dx]=p_dot(p0,psol,'derivatives',[0,1]);
    Jph_dx=p_Jph_dx.profile(:)';
    res=[res;resph];
    nph=1;
    Jstruc.rows.phase=size(J,1)+(1:nph);
    J=[J; Jph_dx(:)',zeros(1,Jstruc.cols.nT),zeros(1,Jstruc.cols.np)];
else
    nph=0;
    Jstruc.rows.phase=size(J,1)+(1:nph);
end
%% add (linear) steplength conditions
nsc=length(options.step_cnd);
Jsc=NaN(nsc,size(J,2));
res_sc=zeros(nsc,1);
for j=1:length(options.step_cnd)
    Jsc(j,Jstruc.cols.xrg)=options.step_cnd(j).profile(:)';
    if Jstruc.cols.nT>0
        Jsc(j,Jstruc.cols.T)=options.step_cnd(j).period;
    end
    if Jstruc.cols.np>0
        Jsc(j,Jstruc.cols.p)=...
            options.step_cnd(j).parameter(free_par);
    end
end
Jstruc.rows.step_cond=size(J,1)+(1:nsc);
J=[J;Jsc];
res=[res;res_sc];
%% add (user-defined) extra conditions
if ~isempty(options.extra_cond)
    ncond=length(options.extra_cond{1});
    Jd=zeros(ncond,size(J,2));
    condi=options.extra_cond{2};
    for k=1:ncond
        Jd(k,1:Jstruc.cols.x(end))=condi(j).profile(:)';
        if Jstruc.cols.nT>0
            Jd(k,Jstruc.cols.T)=condi(k).period;
        end
        if Jstruc.cols.np>0
            Jd(k,Jstruc.np)=condi(k).parameter(free_par);
        end
    end
    res=[res;options.extra_cond{1}];
    Jstruc.rows.extra_cond=size(J,1)+(1:ncond);
    J=[J;Jd];
else
    Jstruc.rows.extra_cond=[];
end    
end
