function df = dir_Holling_rhs_1(xx,par,dxx,dpar,ind)
% 1st directional derivative of HollingTanner rhs
%    05-Mar-2017 17:23:08
%%
x1=xx(1,1,:);
x2=xx(1,2,:);
y1=xx(2,1,:);
y2=xx(2,2,:);
%%
x1_dev=dxx(1,1,:);
x2_dev=dxx(1,2,:);
y1_dev=dxx(2,1,:);
y2_dev=dxx(2,2,:);
%%
a=par(ind.a);
m=par(ind.m);
beta=par(ind.beta);
delta=par(ind.delta);
%%
a_dev=dpar(ind.a);
m_dev=dpar(ind.m);
beta_dev=dpar(ind.beta);
delta_dev=dpar(ind.delta);
h_dev=dpar(ind.h);

t2 = m_dev+x1_dev;
t3 = a.*y1;
t4 = t3+x1;
t5 = 1.0./t4;
df(1,:) = -h_dev-t2.*(m+x1)-t2.*(m+x1-1.0)-t5.*x1.*y1_dev-t5.*x1_dev.*y1+1.0./t4.^2.*x1.*y1.*(x1_dev+a.*y1_dev+a_dev.*y1);
t6 = 1.0./x2;
t7 = beta-t6.*y2;
df(2,:) = delta.*y1.*(beta_dev-t6.*y2_dev+1.0./x2.^2.*x2_dev.*y2)+delta_dev.*t7.*y1+delta.*t7.*y1_dev;
end
