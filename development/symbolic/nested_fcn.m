function fdeqs = nested_fcn(in1,v1,v2,v1_1,v2_1,v1_2,v2_2,v1_3,v2_3,v1_4,v2_4,v1_5,v2_5)
%NESTED_FCN
%    FDEQS = NESTED_FCN(IN1,V1,V2,V1_1,V2_1,V1_2,V2_2,V1_3,V2_3,V1_4,V2_4,V1_5,V2_5)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    03-Mar-2017 16:58:25

xeq1 = in1(1,:);
t2 = v1_1(-xeq1);
t3 = v1(-xeq1);
t4 = t3.^2;
fdeqs = [t4.*v1_2(-xeq1).*-3.0-t2.^2.*t3.*6.0;t4.*v2_2(-xeq1).*-3.0-t2.*t3.*v2_1(-xeq1).*6.0];
