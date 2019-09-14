function [r,J]=m_sys_cond_RWFold(p,orig_cond,dim,ind_rho)
%% constraints used for extended DDE in fold continuation of relative equilibria
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%

%% Orig
p_orig=p; % point strut
p_orig.x=p.x(1:dim); % select actual guess/soln
p_orig.parameter=p_orig.parameter(1:min(ind_rho)-1); % lose rho
[r_orig,J_orig]=orig_cond(p_orig);
nuser=length(r_orig);

numRho = numel(ind_rho);
% This number is used throughout the function. It's equal to
% the number of omegas used in rot_sym conditions.

for i=1:nuser
    J_orig(i).x=[J_orig(i).x;zeros(dim,1)]; % keep orig J, tack on 0s after so you have [x;v]
    J_orig(i).parameter=[J_orig(i).parameter,...
        zeros(1,numRho)]; % put a zero after for rho
end

%% Phas, ORIG from rot_sym
% rphas=0;
% Jphas=J_orig(end);
% Jphas.x=[zeros(dim,1);Jphas.x(1:dim)];

%% I AM GOING TO TRY ADDING TWO phas
% rphas1=0;
% Jphas1=J_orig(end-1); 
% Jphas1.x=[zeros(dim,1);Jphas1.x(1:dim)]; 
% 
% rphas2=0;
% Jphas2=J_orig(end); % Grab the LAST J_orig (related to DDE-BIF rot_sym package.)
% Jphas2.x=[zeros(dim,1);Jphas2.x(1:dim)]; % zeros for x guess, old guess in Jphas for v

%% phas: This is with numel(ind_rho) aka arbitrary amount of omegas
% In this case, numRho ~= 1.
rphasContainer = cell(numRho);
Jphascontainer = cell(numRho);
for i = 1:numRho
    % Create temporary values related to that particular rho.
    rphas = 0;
    Jphas = J_orig(end-numRho+i);
    % This works because if there are two rhos:
    % Jphas(1) = J_orig(end-2+1) aka second to last
    % Jphas(2) = J_orig(end-2+2) aka last
    Jphas.x = [zeros(dim,1);Jphas.x(1:dim)];
    
    % Store temporary values in their containers.
    rphasContainer{i} = rphas;
    Jphascontainer{i} = Jphas;
end

%% Norm
rnorm=sum(p.x(dim+1:end).^2)+sum(p.parameter(ind_rho).^2)-1;
%square each value of v in p, then add each rho^2, then -1.

Jnorm=p_axpy(0,p,[]); % x = 0, par from p
% I assume this is because I am taking each partial derivative of the norm
% squared? I don't know why I would be doing that though...
Jnorm.x(dim+1:end)=2*p.x(dim+1:end); %Give Jnorm, 2*v from p
Jnorm.parameter(ind_rho)=2*p.parameter(ind_rho);
% Since numel(ind_rho) ~= 1: This COULD mean multiple assignments.

%% Return
% From rot_sym
% r=[r_orig(:);rphas;rnorm];
% J=[J_orig(:);Jphas;Jnorm];

% For Two...
% r=[r_orig(:);rphas1;rphas2;rnorm];
% J=[J_orig(:);Jphas1;Jphas2;Jnorm];

% For an arbitrary amount of omegas/numel(ind_rho)
r = r_orig(:);
J = J_orig(:);

% Preallocate, note: numel(r_orig) == numel(J_orig)
r_orig_length = numel(r_orig);
r = [r;zeros(numRho,1)];
J = [J;repmat(struct('kind','stst','parameter',0,'x',0),numRho,1)];
    % repmat is an efficient way to generate structures
% Assign values for return. Overwrite preallocation above.
for i = 1:numRho;
    r(r_orig_length + i) = rphasContainer{i};
    J(r_orig_length + i) = Jphascontainer{i};
end
r=[r;rnorm];
J=[J;Jnorm];

end
