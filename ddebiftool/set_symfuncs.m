function funcs=set_symfuncs(fcnrhs,varargin)
%% convert output of symbolic toolbox into funcs structure recognized by DDE-Biftool
%
%% Inputs
%
% |fcnrhs|: function handle of function created by dde_sym2funcs(...)
%
%% Optional inputs (name value pairs)
%
% * |'sys_tau'|: some inputs are passed on to |set_funcs|. In particular,
% |sys_tau| may need to be specified, if the delays are parameters.
% * |'x_vectorized'| (|true|): if resulting functions can be called with 3d
% array of |xx| arguments [dimension x (number of delays +1) x nvec].
% Typically, the symbolic toolbox generated functionsa permit this.
% * |'p_vectorized'| (|true|): can all the nvec different |xx| inputs also
% have different parameter values? If yes, then the function gets called
% with a parameter array of size 1 x np x nvec. 
%
%% Output
%
% * |funcs|: structure containing system r.h.s, delays, derivatives and
% other information.
%
% $Id: set_symfuncs.m 348 2019-06-19 13:09:01Z jansieber $
%%
default={'x_vectorized',true,'p_vectorized',true}; % change default on vectorization
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if ischar(fcnrhs)
    f=str2func(fcnrhs);
else
    f=fcnrhs;
end
directional_derivative=f('directional_derivative');
if directional_derivative
    sys_rhs=@(x,p)dde_sym_rhs_wrap(f,0,x,p,x,p);
else
    sys_rhs=@(x,p)dde_sym_rhs_wrap(f,0,x,p,{},{});
end
maxorder=f('maxorder');
if directional_derivative
    for order=maxorder:-1:1
        sys_dirderi{order}=@(x,p,dx,dp)dde_sym_rhs_wrap(f,order,x,p,dx,dp);
    end
    sys_mfderi={};
else
    for order=maxorder:-1:1
        sys_dirderi{order}=@(x,p,dx,dp)dde_sym_rhs_wrap(f,order,x,p,dx,dp);
    end
    sys_mfderi=@(xx,p,varargin)wrap_mfderi(f,xx,p,varargin);
end
sys_deri=@(x,p,nx,np,v)dde_gen_deriv(sys_dirderi,x,p,nx,np,v,directional_derivative);
%% construct functions specifying variable delay if requested
sd_delay=f('tp_del');
if sd_delay
    ntau=f('ntau');
    sys_ntau=@()ntau;
    sys_tau=@(itau,x,p)dde_sym_tau_wrap(f,itau,0,x,p,x,p);
    for order=maxorder:-1:1
        sys_dirdtau{order}=@(itau,x,p,dx,dp)dde_sym_tau_wrap(...
        f,itau,order,x,p,dx,dp);
    end
    dfunc=@(itau){...
        @(xa,pa,dx,dp)sys_dirdtau{1}(itau,xa,pa,dx,dp),...
        @(xa,pa,dx,dp)sys_dirdtau{2}(itau,xa,pa,dx,dp)};
    sys_dtau=@(itau,x,p,nx,np)dde_gen_deriv(dfunc(itau),x,p,nx,np,[]);
    tau_args={'sys_tau',sys_tau,'sys_dtau',sys_dtau,'sys_ntau',sys_ntau,...
        'sys_dirdtau',sys_dirdtau};
else
    tau_args={};
end
funcs=set_funcs('sys_rhs',sys_rhs,'sys_deri',sys_deri,'sys_dirderi',sys_dirderi,...
    tau_args{:},'sys_mfderi',sys_mfderi,...
    'x_vectorized',options.x_vectorized,'p_vectorized',options.p_vectorized,...
    pass_on{:});
end
%%
function y=wrap_mfderi(f,xx,p,devs)
order=length(devs);
np=length(p);
[n,ntau]=size(xx);
dx=zeros(n,ntau,order);
dp=zeros(np,order);
for i=1:order
   devi=reshape(devs{i},n+np,ntau);
   dx(:,:,i)=devi(1:n,:);
   dp(:,i)=devi(n+1:end,1);
end
xx=xx(1:n,:);
y=dde_sym_rhs_wrap(f,order,xx,p,dx,dp);
end
