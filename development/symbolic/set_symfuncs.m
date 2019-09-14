function funcs=set_symfuncs(fcnrhs,varargin)
%% convert output of symbolic toolbox into funcs structure recognized by DDE-Biftool
%
% $Id: set_symfuncs.m 174 2017-03-10 21:59:17Z jansieber $
%%
default={'x_vectorized',true}; % change default on vectorization
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if ischar(fcnrhs)
    f=str2func(fcnrhs);
else
    f=fcnrhs;
end
sys_rhs=@(x,p)dde_sym_rhs_wrap(f,0,x,p,x,p);
sys_dirderi=@(order,x,p,dx,dp)dde_sym_rhs_wrap(f,order,x,p,dx,dp);
sys_deri=@(x,p,nx,np,v)dde_gen_deriv(sys_dirderi,x,p,nx,np,v);
%% construct functions specifying variable delay if requested
sd_delay=f('tp_del');
if sd_delay
    mf_dxlength=f('mf_dxlength');
    ntau=f('ntau');
    sys_ntau=@()ntau;
    sys_tau=@(itau,x,p)dde_sym_tau_wrap(f,itau,0,x,p,x,p);
    sys_dirdtau=@(itau,order,x,p,dx,dp)dde_sym_tau_wrap(...
        f,itau,order,x,p,dx,dp);
    sys_dtau=@(itau,x,p,nx,np)dde_gen_deriv(...
        @(order,xa,pa,dxa,dpa)sys_dirdtau(itau,order,xa,pa,dxa,dpa),x,p,nx,np,[]);
    sys_dirmf=@(order,x,v,p)dde_sym_sd_mfderi_wrap(f,order,mf_dxlength,x,p,v);
    tau_args={'sys_tau',sys_tau,'sys_dtau',sys_dtau,'sys_ntau',sys_ntau,...
        'sys_dirdtau',sys_dirdtau,'sys_dirmf',sys_dirmf};
else
    tau_args={};
end
funcs=set_funcs('sys_rhs',sys_rhs,'sys_deri',sys_deri,'sys_dirderi',sys_dirderi,...
    tau_args{:},'x_vectorized',options.x_vectorized,pass_on{:});
end
