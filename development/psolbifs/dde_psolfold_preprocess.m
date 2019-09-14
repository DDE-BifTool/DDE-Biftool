function data=dde_psolfold_preprocess(data)
data.orig=data;
%% convert points
data.step_cnd=arrayfun(@dde_coll_convert,data.orig.step_cnd);
pt=data.orig.point;
[data.point,data.conversion]=dde_coll_convert(data.orig.point);
data.previous=dde_coll_convert(data.orig.previous);
data.method.extra_condition=true;
npar=length(pt.parameter);
npar_var=length(pt.parameter_var);
if ~data.funcs.tp_del
    npar_del=length(pt.parameter_delay);
    ipdelays=npar+npar_var+(1:npar_del);
    itau=data.orig.funcs.sys_tau();
    ntau1=length(itau)+1;
    sys_tau=@()[itau(:)',ipdelays];
    [dum,delay_free]=ismember(itau,data.orig.free_par); %#ok<ASGLU>
    delay_free=delay_free(delay_free>0);
    data.free_par=[data.orig.free_par,npar+(1:npar_var),ipdelays];
    dtaudp=sparse((2:ntau1)',itau(:),ones(ntau1-1,1),ntau1,npar);
    dtaudp=full(dtaudp);
    sys_rhs=@(x,p)sys_rhs_psolfold(x,p(1:npar),p(npar+(1:npar_var)),[],[],[],...
        data.orig.funcs,pt.free_par,data.orig.free_par,dtaudp);
    if ~isfield(data.orig.method,'hdev') || data.orig.method.hdev==0
        sys_dirderi={'sys_dirderi',{@(x,p,dx,dp)sys_rhs_psolfold(...
            x,p(1:npar),p(npar+(1:npar_var)),dx,dp(1:npar),dp(npar+(1:npar_var)),...
            data.orig.funcs,pt.free_par,data.orig.free_par,dtaudp)}};
    else
        sys_dirderi={'hjac',@(ord)data.orig.method.hdev.^(1./ord)};
    end
    sys_cond=@(p)sys_cond_psolfold(p,...
        data.orig.funcs,data.conversion,data.method);
    data.funcs=set_funcs(...
        'sys_rhs',sys_rhs,...
        sys_dirderi{:},...
        'sys_tau',sys_tau,...
        'sys_cond',sys_cond,...
        'x_vectorized',true);
    data.funcs.delayed_derivs=[zeros(1,length(itau)),...
        ones(1,length(itau)+1)];
else
    %% state-dependent delay
    ntau1=data.orig.sys_ntau()+1;
    tau=data.orig.funcs.sys_tau;
    dtau=data.orig.funcs.sys_dtau;    
    sys_tau=@(k,x,p)sys_tau_psolbif(k,ntau1,tau,x,p);
    sys_dtau=@(k,x,p,nx,np)sys_tau_psolbif(k,ntau1,dtau,x,p,nx,np);
    sys_rhs=@(x,p)sys_rhs_psolfold(x,p(1:npar),p(npar+(1:npar_var)),[],[],[],...
        data.orig.funcs,pt.free_par,data.orig.free_par,[]);
end
end
