function varargout=sym_poscontrol(action,varargin)
%% Automatically generated with matlabFunction
% 
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'ntau'
   varargout{1}=3;
   return
  case 'tp_del'
   varargout{1}=1;
   return
  case 'maxorder'
   varargout{1}=5;
   return
end
ind=varargin{1};
order=varargin{2};
nout=varargin{3};
f=str2func(sprintf('sym_poscontrol_%s_%d_%d',action,ind,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{4:end});




function [out1,out2] = sym_poscontrol_rhs_1_0(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_RHS_1_0
%    [OUT1,OUT2] = SYM_POSCONTROL_RHS_1_0(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:30:23

out1 = K.*c.*(s0-s2).*(1.0./2.0);
if nargout > 1
    out2 = (x1+x3-c.*s1-K.*c.*(s0.*-2.0+s2+s4).*(1.0./2.0))./(c-K.*c.*(s0-s4).*(1.0./2.0));
end


function [out1,out2] = sym_poscontrol_rhs_1_1(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_RHS_1_1
%    [OUT1,OUT2] = SYM_POSCONTROL_RHS_1_1(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:30:43

t2 = s0-s2;
out1 = K.*c_dev.*t2.*(1.0./2.0)+K_dev.*c.*t2.*(1.0./2.0)+K.*c.*(s0_dev-s2_dev).*(1.0./2.0);
if nargout > 1
    t6 = s0.*2.0;
    t3 = s2+s4-t6;
    t4 = s0-s4;
    t5 = c-K.*c.*t4.*(1.0./2.0);
    out2 = -(-x1_dev-x3_dev+c.*s1_dev+c_dev.*s1+K.*c.*(s0_dev.*-2.0+s2_dev+s4_dev).*(1.0./2.0)+K.*c_dev.*t3.*(1.0./2.0)+K_dev.*c.*t3.*(1.0./2.0))./t5+1.0./t5.^2.*(x1+x3-c.*s1-K.*c.*t3.*(1.0./2.0)).*(-c_dev+K.*c_dev.*t4.*(1.0./2.0)+K_dev.*c.*t4.*(1.0./2.0)+K.*c.*(s0_dev-s4_dev).*(1.0./2.0));
end


function [out1,out2] = sym_poscontrol_rhs_1_2(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_RHS_1_2
%    [OUT1,OUT2] = SYM_POSCONTROL_RHS_1_2(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:30:48

t2 = s0_dev-s2_dev;
out1 = K.*c_dev.*t2+K_dev.*c.*t2+K_dev.*c_dev.*(s0-s2);
if nargout > 1
    t10 = s0_dev.*2.0;
    t3 = s2_dev+s4_dev-t10;
    t4 = s0-s4;
    t8 = K.*c.*t4.*(1.0./2.0);
    t5 = c-t8;
    t6 = s0_dev-s4_dev;
    t11 = s0.*2.0;
    t7 = s2+s4-t11;
    t9 = 1.0./t5.^2;
    t12 = K.*c.*t6.*(1.0./2.0);
    t13 = K.*c_dev.*t4.*(1.0./2.0);
    t14 = K_dev.*c.*t4.*(1.0./2.0);
    t15 = -c_dev+t12+t13+t14;
    t16 = x1+x3-c.*s1-K.*c.*t7.*(1.0./2.0);
    out2 = -(c_dev.*s1_dev.*2.0+K.*c_dev.*t3+K_dev.*c.*t3+K_dev.*c_dev.*t7)./t5+t9.*t16.*(K.*c_dev.*t6+K_dev.*c.*t6+K_dev.*c_dev.*t4)-t9.*t15.*(-x1_dev-x3_dev+c.*s1_dev+c_dev.*s1+K.*c.*t3.*(1.0./2.0)+K.*c_dev.*t7.*(1.0./2.0)+K_dev.*c.*t7.*(1.0./2.0)).*2.0+1.0./t5.^3.*t15.^2.*t16.*2.0;
end


function [out1,out2] = sym_poscontrol_rhs_1_3(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_RHS_1_3
%    [OUT1,OUT2] = SYM_POSCONTROL_RHS_1_3(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:30:49

out1 = K_dev.*c_dev.*(s0_dev-s2_dev).*3.0;
if nargout > 1
    t14 = s0_dev.*2.0;
    t2 = s2_dev+s4_dev-t14;
    t3 = s0-s4;
    t12 = K.*c.*t3.*(1.0./2.0);
    t4 = c-t12;
    t5 = s0_dev-s4_dev;
    t6 = K.*c.*t5.*(1.0./2.0);
    t7 = K.*c_dev.*t3.*(1.0./2.0);
    t8 = K_dev.*c.*t3.*(1.0./2.0);
    t9 = -c_dev+t6+t7+t8;
    t10 = t9.^2;
    t15 = s0.*2.0;
    t11 = s2+s4-t15;
    t13 = 1.0./t4.^2;
    t16 = c.*s1_dev;
    t17 = c_dev.*s1;
    t18 = K.*c.*t2.*(1.0./2.0);
    t19 = K.*c_dev.*t11.*(1.0./2.0);
    t20 = K_dev.*c.*t11.*(1.0./2.0);
    t21 = t16+t17+t18+t19+t20-x1_dev-x3_dev;
    t22 = 1.0./t4.^3;
    t23 = K.*c_dev.*t5;
    t24 = K_dev.*c.*t5;
    t25 = K_dev.*c_dev.*t3;
    t26 = t23+t24+t25;
    t28 = c.*s1;
    t29 = K.*c.*t11.*(1.0./2.0);
    t27 = -t28-t29+x1+x3;
    out2 = t9.*t13.*(c_dev.*s1_dev.*2.0+K.*c_dev.*t2+K_dev.*c.*t2+K_dev.*c_dev.*t11).*-3.0-t10.*t21.*t22.*6.0-t13.*t21.*t26.*3.0+t9.*t22.*t26.*t27.*6.0-(K_dev.*c_dev.*t2.*3.0)./t4+1.0./t4.^4.*t9.*t10.*t27.*6.0+K_dev.*c_dev.*t5.*t13.*t27.*3.0;
end


function [out1,out2] = sym_poscontrol_rhs_1_4(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_RHS_1_4
%    [OUT1,OUT2] = SYM_POSCONTROL_RHS_1_4(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:30:51

out1 = 0.0;
if nargout > 1
    t2 = s0-s4;
    t6 = s0_dev-s4_dev;
    t16 = K.*c.*t6.*(1.0./2.0);
    t17 = K.*c_dev.*t2.*(1.0./2.0);
    t18 = K_dev.*c.*t2.*(1.0./2.0);
    t3 = -c_dev+t16+t17+t18;
    t4 = t3.^2;
    t9 = K.*c.*t2.*(1.0./2.0);
    t5 = c-t9;
    t19 = s0_dev.*2.0;
    t7 = s2_dev+s4_dev-t19;
    t14 = s0.*2.0;
    t8 = s2+s4-t14;
    t10 = K.*c_dev.*t6;
    t11 = K_dev.*c.*t6;
    t12 = K_dev.*c_dev.*t2;
    t13 = t10+t11+t12;
    t33 = c.*s1;
    t34 = K.*c.*t8.*(1.0./2.0);
    t15 = -t33-t34+x1+x3;
    t20 = 1.0./t5.^3;
    t21 = c_dev.*s1_dev.*2.0;
    t22 = K.*c_dev.*t7;
    t23 = K_dev.*c.*t7;
    t24 = K_dev.*c_dev.*t8;
    t25 = t21+t22+t23+t24;
    t26 = c.*s1_dev;
    t27 = c_dev.*s1;
    t28 = K.*c.*t7.*(1.0./2.0);
    t29 = K.*c_dev.*t8.*(1.0./2.0);
    t30 = K_dev.*c.*t8.*(1.0./2.0);
    t31 = t26+t27+t28+t29+t30-x1_dev-x3_dev;
    t32 = 1.0./t5.^4;
    t35 = 1.0./t5.^2;
    out2 = t13.^2.*t15.*t20.*6.0+t4.^2.*1.0./t5.^5.*t15.*2.4e1-t4.*t20.*t25.*1.2e1-t13.*t25.*t35.*6.0+t4.*t13.*t15.*t32.*3.6e1-t3.*t13.*t20.*t31.*2.4e1-t3.*t4.*t31.*t32.*2.4e1-K_dev.*c_dev.*t3.*t7.*t35.*1.2e1-K_dev.*c_dev.*t6.*t31.*t35.*1.2e1+K_dev.*c_dev.*t3.*t6.*t15.*t20.*2.4e1;
end


function [out1,out2] = sym_poscontrol_rhs_1_5(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_RHS_1_5
%    [OUT1,OUT2] = SYM_POSCONTROL_RHS_1_5(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:30:53

out1 = 0.0;
if nargout > 1
    t2 = s0-s4;
    t7 = s0_dev-s4_dev;
    t8 = K.*c.*t7.*(1.0./2.0);
    t9 = K.*c_dev.*t2.*(1.0./2.0);
    t10 = K_dev.*c.*t2.*(1.0./2.0);
    t3 = -c_dev+t8+t9+t10;
    t4 = t3.^2;
    t5 = t4.^2;
    t13 = K.*c.*t2.*(1.0./2.0);
    t6 = c-t13;
    t12 = s0.*2.0;
    t11 = s2+s4-t12;
    t15 = s0_dev.*2.0;
    t14 = s2_dev+s4_dev-t15;
    t24 = K.*c_dev.*t7;
    t25 = K_dev.*c.*t7;
    t26 = K_dev.*c_dev.*t2;
    t16 = t24+t25+t26;
    t17 = c.*s1_dev;
    t18 = c_dev.*s1;
    t19 = K.*c.*t14.*(1.0./2.0);
    t20 = K.*c_dev.*t11.*(1.0./2.0);
    t21 = K_dev.*c.*t11.*(1.0./2.0);
    t22 = t17+t18+t19+t20+t21-x1_dev-x3_dev;
    t23 = 1.0./t6.^4;
    t27 = 1.0./t6.^3;
    t28 = c_dev.*s1_dev.*2.0;
    t29 = K.*c_dev.*t14;
    t30 = K_dev.*c.*t14;
    t31 = K_dev.*c_dev.*t11;
    t32 = t28+t29+t30+t31;
    t33 = t16.^2;
    t36 = c.*s1;
    t37 = K.*c.*t11.*(1.0./2.0);
    t34 = -t36-t37+x1+x3;
    t35 = 1.0./t6.^5;
    t38 = 1.0./t6.^2;
    out2 = t5.*t22.*t35.*-1.2e2-t22.*t27.*t33.*3.0e1-t3.*t4.*t23.*t32.*6.0e1-t4.*t16.*t22.*t23.*1.8e2-t3.*t16.*t27.*t32.*6.0e1+t3.*t23.*t33.*t34.*9.0e1+t3.*t5.*1.0./t6.^6.*t34.*1.2e2-K_dev.*c_dev.*t4.*t14.*t27.*6.0e1-K_dev.*c_dev.*t14.*t16.*t38.*3.0e1-K_dev.*c_dev.*t7.*t32.*t38.*3.0e1+t3.*t4.*t16.*t34.*t35.*2.4e2-K_dev.*c_dev.*t3.*t7.*t22.*t27.*1.2e2+K_dev.*c_dev.*t4.*t7.*t23.*t34.*1.8e2+K_dev.*c_dev.*t7.*t16.*t27.*t34.*6.0e1;
end


function out1 = sym_poscontrol_tau_1_0(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_1_0
%    OUT1 = SYM_POSCONTROL_TAU_1_0(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:06

out1 = tau0;


function out1 = sym_poscontrol_tau_2_0(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_2_0
%    OUT1 = SYM_POSCONTROL_TAU_2_0(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:08

out1 = s1;


function out1 = sym_poscontrol_tau_3_0(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_3_0
%    OUT1 = SYM_POSCONTROL_TAU_3_0(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:08

out1 = s1+tau0;


function out1 = sym_poscontrol_tau_1_1(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_1_1
%    OUT1 = SYM_POSCONTROL_TAU_1_1(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:09

out1 = tau0_dev;


function out1 = sym_poscontrol_tau_2_1(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_2_1
%    OUT1 = SYM_POSCONTROL_TAU_2_1(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:09

out1 = s1_dev;


function out1 = sym_poscontrol_tau_3_1(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_3_1
%    OUT1 = SYM_POSCONTROL_TAU_3_1(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:10

out1 = s1_dev+tau0_dev;


function out1 = sym_poscontrol_tau_1_2(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_1_2
%    OUT1 = SYM_POSCONTROL_TAU_1_2(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:10

out1 = 0.0;


function out1 = sym_poscontrol_tau_2_2(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_2_2
%    OUT1 = SYM_POSCONTROL_TAU_2_2(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:10

out1 = 0.0;


function out1 = sym_poscontrol_tau_3_2(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_3_2
%    OUT1 = SYM_POSCONTROL_TAU_3_2(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:11

out1 = 0.0;


function out1 = sym_poscontrol_tau_1_3(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_1_3
%    OUT1 = SYM_POSCONTROL_TAU_1_3(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:11

out1 = 0.0;


function out1 = sym_poscontrol_tau_2_3(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_2_3
%    OUT1 = SYM_POSCONTROL_TAU_2_3(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:11

out1 = 0.0;


function out1 = sym_poscontrol_tau_3_3(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_3_3
%    OUT1 = SYM_POSCONTROL_TAU_3_3(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:12

out1 = 0.0;


function out1 = sym_poscontrol_tau_1_4(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_1_4
%    OUT1 = SYM_POSCONTROL_TAU_1_4(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:12

out1 = 0.0;


function out1 = sym_poscontrol_tau_2_4(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_2_4
%    OUT1 = SYM_POSCONTROL_TAU_2_4(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:13

out1 = 0.0;


function out1 = sym_poscontrol_tau_3_4(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_3_4
%    OUT1 = SYM_POSCONTROL_TAU_3_4(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:13

out1 = 0.0;


function out1 = sym_poscontrol_tau_1_5(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_1_5
%    OUT1 = SYM_POSCONTROL_TAU_1_5(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:13

out1 = 0.0;


function out1 = sym_poscontrol_tau_2_5(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_2_5
%    OUT1 = SYM_POSCONTROL_TAU_2_5(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:13

out1 = 0.0;


function out1 = sym_poscontrol_tau_3_5(x1,s1,x2,s2,x3,s3,x4,s4,tau0,s0,K,c,x1_dev,s1_dev,x2_dev,s2_dev,x3_dev,s3_dev,x4_dev,s4_dev,tau0_dev,s0_dev,K_dev,c_dev)
%SYM_POSCONTROL_TAU_3_5
%    OUT1 = SYM_POSCONTROL_TAU_3_5(X1,S1,X2,S2,X3,S3,X4,S4,TAU0,S0,K,C,X1_DEV,S1_DEV,X2_DEV,S2_DEV,X3_DEV,S3_DEV,X4_DEV,S4_DEV,TAU0_DEV,S0_DEV,K_DEV,C_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    13-Jul-2017 21:31:14

out1 = 0.0;

