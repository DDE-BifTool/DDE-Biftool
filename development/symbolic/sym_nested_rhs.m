function varargout=sym_nested_rhs(ind,order,nout,varargin)
%% Automatically generated with matlabFunction
% 
%#ok<*DEFNU,*INUSD,*INUSL>
f=str2func(sprintf('sym_nested_rhs_%d_%d',ind,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{:});


function [out_1,out_2] = sym_nested_rhs_1_0(xx1_1,xx2_1,xx1_2,xx2_2,p1,p2,xx1_1_dev,xx2_1_dev,xx1_2_dev,xx2_2_dev,p1_dev,p2_dev)
%SYM_NESTED_RHS_1_0
%    [OUT_1,OUT_2] = SYM_NESTED_RHS_1_0(XX1_1,XX2_1,XX1_2,XX2_2,P1,P2,XX1_1_DEV,XX2_1_DEV,XX1_2_DEV,XX2_2_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    05-Mar-2017 00:01:56

out_1 = p1-xx1_2;
if nargout > 1
    out_2 = -xx2_2;
end


function [out_1,out_2] = sym_nested_rhs_1_1(xx1_1,xx2_1,xx1_2,xx2_2,p1,p2,xx1_1_dev,xx2_1_dev,xx1_2_dev,xx2_2_dev,p1_dev,p2_dev)
%SYM_NESTED_RHS_1_1
%    [OUT_1,OUT_2] = SYM_NESTED_RHS_1_1(XX1_1,XX2_1,XX1_2,XX2_2,P1,P2,XX1_1_DEV,XX2_1_DEV,XX1_2_DEV,XX2_2_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    05-Mar-2017 00:01:57

out_1 = p1_dev-xx1_2_dev;
if nargout > 1
    out_2 = -xx2_2_dev;
end


function [out_1,out_2] = sym_nested_rhs_1_2(xx1_1,xx2_1,xx1_2,xx2_2,p1,p2,xx1_1_dev,xx2_1_dev,xx1_2_dev,xx2_2_dev,p1_dev,p2_dev)
%SYM_NESTED_RHS_1_2
%    [OUT_1,OUT_2] = SYM_NESTED_RHS_1_2(XX1_1,XX2_1,XX1_2,XX2_2,P1,P2,XX1_1_DEV,XX2_1_DEV,XX1_2_DEV,XX2_2_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    05-Mar-2017 00:01:58

out_1 = 0.0;
if nargout > 1
    out_2 = 0.0;
end


function [out_1,out_2] = sym_nested_rhs_1_3(xx1_1,xx2_1,xx1_2,xx2_2,p1,p2,xx1_1_dev,xx2_1_dev,xx1_2_dev,xx2_2_dev,p1_dev,p2_dev)
%SYM_NESTED_RHS_1_3
%    [OUT_1,OUT_2] = SYM_NESTED_RHS_1_3(XX1_1,XX2_1,XX1_2,XX2_2,P1,P2,XX1_1_DEV,XX2_1_DEV,XX1_2_DEV,XX2_2_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    05-Mar-2017 00:01:58

out_1 = 0.0;
if nargout > 1
    out_2 = 0.0;
end


function [out_1,out_2] = sym_nested_rhs_1_4(xx1_1,xx2_1,xx1_2,xx2_2,p1,p2,xx1_1_dev,xx2_1_dev,xx1_2_dev,xx2_2_dev,p1_dev,p2_dev)
%SYM_NESTED_RHS_1_4
%    [OUT_1,OUT_2] = SYM_NESTED_RHS_1_4(XX1_1,XX2_1,XX1_2,XX2_2,P1,P2,XX1_1_DEV,XX2_1_DEV,XX1_2_DEV,XX2_2_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    05-Mar-2017 00:01:58

out_1 = 0.0;
if nargout > 1
    out_2 = 0.0;
end


function [out_1,out_2] = sym_nested_rhs_1_5(xx1_1,xx2_1,xx1_2,xx2_2,p1,p2,xx1_1_dev,xx2_1_dev,xx1_2_dev,xx2_2_dev,p1_dev,p2_dev)
%SYM_NESTED_RHS_1_5
%    [OUT_1,OUT_2] = SYM_NESTED_RHS_1_5(XX1_1,XX2_1,XX1_2,XX2_2,P1,P2,XX1_1_DEV,XX2_1_DEV,XX1_2_DEV,XX2_2_DEV,P1_DEV,P2_DEV)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    05-Mar-2017 00:01:58

out_1 = 0.0;
if nargout > 1
    out_2 = 0.0;
end
