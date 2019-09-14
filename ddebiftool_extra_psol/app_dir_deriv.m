function dy=app_dir_deriv(f,x0,deviation,ind,h)
%% 2nd order approx directional derivative in ith column of x0
%
% $Id: app_dir_deriv.m 309 2018-10-28 19:02:42Z jansieber $
%
xd=x0;
xd(:,ind,:)=x0(:,ind,:)+h*deviation;
dy1=f(xd);
xd(:,ind,:)=x0(:,ind,:)-h*deviation;
dy2=f(xd);
dy=(dy1-dy2)/(2*h);
end