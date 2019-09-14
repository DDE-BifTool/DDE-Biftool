function y=m_sys_rhs_RWFold(x,p,rho,orig_rhs,dim,hbif)
%% rhs of extended DDE for fold of relative equilibria
%
% x extended state (orig state and null vector)
% p user parameters
% rho rotation frequency, double = old behavior, cell = new behavior
% orig_rhs user r.h.s
% dim original system dimension
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
%% Now adding another omega!!
x0=x(1:dim,:); % the actual guess/soln vector
v=x(dim+1:end,:); % the nullvector
y0=orig_rhs(x0,p);
dfx=@(x0,dev,ind)app_dir_deriv(@(x)orig_rhs(x,p),x0,dev,ind,hbif);
%      2nd order approx directional derivative in ith column of x0
% 
%      dy = app_dir_deriv(f,x0,deviation,ind,h)
%      xd=x0;
%      xd(:,ind,:)=x0(:,ind,:)+h*deviation;
%      dy1=f(xd);
%      xd(:,ind,:)=x0(:,ind,:)-h*deviation;
%      dy2=f(xd);
%      dy=(dy1-dy2)/(2*h);

%% add partial derivatives of all terms
y1=0*y0; % Make a vector that is the right size but all zeros.
for i=1:size(x0,2)
    y1=y1+dfx(x0,v(:,i),i);
end

% HOW IT WORKS WITH 2 OMEGA, WITH THE LAST TWO PARAMS (... ,omega1, omega2)
%
%     % For omega 1
%     dfom1=app_dir_deriv(...
%         @(p)orig_rhs(x0,p),... % function to provide, devations are in par
%         p,...                  % parameter set
%         1,...                  % deviation = 1 * hbif
%         length(p)-1,...        % select the parameter which is perturbed, omega
%         hbif);                 % this is the size of the deviation.
% 
%     % For omega 2
%     dfom2=app_dir_deriv(...
%         @(p)orig_rhs(x0,p),... % function to provide, devations are in par
%         p,...                  % parameter set
%         1,...                  % deviation = 1 * hbif
%         length(p),...          % select the parameter which is perturbed, omega
%         hbif);                 % this is the size of the deviation.

numRho = numel(rho);
dfom = cell(numel(rho),1);
% Create a dir deriv for each rho parameter we have. Organize them
% into a cell, so they can multiple by the parameter later on.
for i = 1:numel(rho)
    TEMPdfom = app_dir_deriv(...
    @(p)orig_rhs(x0,p),... % function to provide, devations are in par
    p,...                  % parameter set
    1,...                  % deviation = 1 * hbif
    length(p)-numRho+i,... % select the parameter which is perturbed, omega
    hbif);                 % this is the size of the deviation.

    dfom{i} = TEMPdfom;
end

% add each rho related deriv
for i = 1:numel(rho)
    y1 = y1 + dfom{i}*rho(i);
end
y=[y0;y1];


end
