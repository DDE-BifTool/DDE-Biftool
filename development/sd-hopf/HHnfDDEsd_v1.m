% In this script file we compute the normal form of the st dep delay diff
% eq:
%
%   u'(t) = -gamma*u(t) - K1*u(t-a1-c*u(t)) - K2*u(t-a2-c*u(t))
%
% on the four dimensional center manifold at a double hopf bifurcation. 
% 
% We follow the derivation reduction on the center manifold and derivation
% of the normal form from the following references.
%
% [1] Shangjiang Guo and Jianhong Wu, "Bifurcation Theory
%     of Functional Differential Equations", 2013
% [2] Y. A. Kuznetsov. Elements of applied bifurcation theory, volume 112 
% of Applied Mathematical Sciences. Springer-Verlag, New York, third 
% edition, 2004.
%
% Our nonlinearity is given by:
%
% N(phi) = K1*(phi(-a1) - phi(-a1 -c phi(0)))
%         + K2*(phi(-a2) - phi(-a2 -c phi(0)))
%
% The flow on the centre manifold is give by:
% 
% z_j' = i*omega_j z_j(t) + conj(D_j)* N(u_t)
%

clearvars;

tic;

syms x y w v;
    
syms omega1 omega2 theta real;

I = sqrt(-1); 

%% Declare Hopf-Hopf point HHj to find:
HHj=1;

%% Compute the expansion of the nonlinearity.

% The basis phi is given by psi = (q_1, conj(q_1), q_2, conj(q_2))
% where 
%
% q_1 = exp(I*omega1*theta)
% q_2 = exp(I*omega2*theta)

phiz = exp(I*omega1*theta)*x + exp(-I*omega1*theta)*y ...
    + exp(I*omega2*theta)*w + exp(-I*omega2*theta)*v;

syms tau gamma K1 K2 c a1 a2 real;

% We compute the expansion of the nonlinearity up to order 3 from the 
% following definition of f(phi*z).

fphiz = K1* ( subs(phiz,theta, -a1) - subs(phiz,theta, -a1 - c* subs(phiz, theta,0)))...
  + K2* ( subs(phiz,theta, -a2) - subs(phiz,theta,  -a2 - c* subs(phiz, theta,0)));

%% Construct the bilinear and trilinear forms with the expansion of the nonlinearity
        
% We will compute the expansion and store the coefficients in the
% matrices F F1 F2 F3 F4. The notation for the coefficients will be that if
% x = xi_1, y = xi_2, w = xi_3, w = xi_4, then
% F(i,j) is the coefficient of xi_i * xi_j, and Fk(i,j) is the coefficient
% of xi_k*xi_i*xi_j for k,i,j in {1,2,3,4}. Notice that this covers all
% the terms of the expansion of the nonlinearity up to order 3.
        
syms F F1 F2 F3 F4
        
% Coefficients of quadratic terms. Stored in matrix F
        
F(1,1) = subs(diff(diff(fphiz,'x'),'x'), {x,y,w,v}, [0,0,0,0]);
F(1,2) = subs(diff(diff(fphiz,'x'),'y'), {x,y,w,v}, [0,0,0,0]);
F(1,3) = subs(diff(diff(fphiz,'x'),'w'), {x,y,w,v}, [0,0,0,0]);
F(1,4) = subs(diff(diff(fphiz,'x'),'v'), {x,y,w,v}, [0,0,0,0]);

F(2,1) = subs(diff(diff(fphiz,'y'),'x'), {x,y,w,v}, [0,0,0,0]);
F(2,2) = subs(diff(diff(fphiz,'y'),'y'), {x,y,w,v}, [0,0,0,0]);
F(2,3) = subs(diff(diff(fphiz,'y'),'w'), {x,y,w,v}, [0,0,0,0]);
F(2,4) = subs(diff(diff(fphiz,'y'),'v'), {x,y,w,v}, [0,0,0,0]);

F(3,1) = subs(diff(diff(fphiz,'w'),'x'), {x,y,w,v}, [0,0,0,0]);
F(3,2) = subs(diff(diff(fphiz,'w'),'y'), {x,y,w,v}, [0,0,0,0]);
F(3,3) = subs(diff(diff(fphiz,'w'),'w'), {x,y,w,v}, [0,0,0,0]);
F(3,4) = subs(diff(diff(fphiz,'w'),'v'), {x,y,w,v}, [0,0,0,0]);

F(4,1) = subs(diff(diff(fphiz,'v'),'x'), {x,y,w,v}, [0,0,0,0]);
F(4,2) = subs(diff(diff(fphiz,'v'),'y'), {x,y,w,v}, [0,0,0,0]);
F(4,3) = subs(diff(diff(fphiz,'v'),'w'), {x,y,w,v}, [0,0,0,0]);
F(4,4) = subs(diff(diff(fphiz,'v'),'v'), {x,y,w,v}, [0,0,0,0]);

%Coefficients of cubic terms. Stored in matrices F1, F2, F3, F4

F1(1,2) = subs(diff(diff(diff(fphiz,'x'),'x'),'y'), {x,y,w,v}, [0,0,0,0]);

F1(2,3) = subs(diff(diff(diff(fphiz,'x'),'y'),'w'), {x,y,w,v}, [0,0,0,0]);

F1(3,4) = subs(diff(diff(diff(fphiz,'x'),'w'),'v'), {x,y,w,v}, [0,0,0,0]);

F3(3,4) = subs(diff(diff(diff(fphiz,'w'),'w'),'v'), {x,y,w,v}, [0,0,0,0]);

%% Construct the nonlinearity in the centre manifold variables

% We compute the vector psi(0), where psi is the basis for the adjoint 
% problem. With
% 
%  psi=(p_1, conj(p_1), p_2, conj(p_2))
%
%   p_1 =D1*exp(I*theta*omega1);
%   p_2 =D1*exp(I*theta*omega2);

D1 = 1/(1 - K1*a1*exp(I*a1*omega1) - K2*a2*exp(I*a2*omega1));

D2 = 1/(1 - K1*a1*exp(I*a1*omega2) - K2*a2*exp(I*a2*omega2));

p_1 =D1*exp(I*theta*omega1);
p_2 =D2*exp(I*theta*omega2);


%% Compute values of the parameters at the Hopf-Hopf point HHj.

% Values for the Hopf-Hopf point.

a1v = 1.3;
a2v = 6;
gammav = 4.75;

ttoc=toc;
[K1v,K2v, omega1v,omega2v] = findHH(HHj);  
tttoc=toc;
cv = 1;

%% Build the bilinear and trilinear forms needed for the normal form computation.

% The formulas are contained in [1]. We try to follow the notation in [1] as 
% closely as we are able without making the notation unpractical.

syms g1_2000 g2_2000 g1_0020 g2_0020 g1_1100 g2_1100 g1_1010 
syms g2_1010 g1_1001 g2_1001 g1_0200 g2_0200 g1_0002 g2_0002 
syms g1_0011 g2_0011 g1_0101 g2_0101 g1_0110 g2_0110 

g1_2000 = conj(D1)*F(1,1);
g2_2000 = conj(D2)*F(1,1);

g1_0020 = conj(D1)*F(3,3); 
g2_0020 = conj(D2)*F(3,3);

g1_1100 = conj(D1)*F(1,2); 
g2_1100 = conj(D2)*F(1,2);

g1_1010 = conj(D1)*F(1,3); 
g2_1010 = conj(D2)*F(1,3);

g1_1001 = conj(D1)*F(1,4); 
g2_1001 = conj(D2)*F(1,4);

g1_0200 = conj(D1)*F(2,2); 
g2_0200 = conj(D2)*F(2,2);

g1_0002 = conj(D1)*F(4,4); 
g2_0002 = conj(D2)*F(4,4);

g1_0011 = conj(D1)*F(3,4); 
g2_0011 = conj(D2)*F(3,4);

g1_0101 = conj(D1)*F(2,4); 
g2_0101 = conj(D2)*F(2,4);

g1_0110 = conj(D1)*F(3,2); 
g2_0110 = conj(D2)*F(3,2);

q_1 = exp(I*omega1*theta);
q_2 = exp(I*omega2*theta);

syms Delta lambda

Delta = lambda + gamma +K1*exp(-a1*lambda) +K2*exp(-a2*lambda);

E1100 = F(1,2)/(subs(Delta, lambda, 0));

E2000 = F(1,1)/(subs(Delta, lambda, 2*I*omega1));

E1010 = F(1,3)/(subs(Delta, lambda, I*(omega1+omega2)));

E1001 = F(1,4)/(subs(Delta, lambda, I*(omega1-omega2)));

E0020 = F(3,3)/(subs(Delta, lambda, 2*I*omega2));

E0011 = F(3,4)/(subs(Delta, lambda, 0));

syms w2000 w1100 w1010 w1001 w0020 w0011 
syms theta real

w2000 = -g1_2000 *exp(I*omega1*theta)/(I*omega1)...
    -conj(g1_0200*exp(I*omega1*theta))/(3*I*omega1)... 
    +g2_2000*exp(I*omega2*theta)/(I*(omega2-2*omega1))...
    -conj(g2_0200*exp(I*omega2*theta))/(I*(omega2+2*omega1)) ...
    +E2000*exp(2*I*omega1*theta);

w1100 = g1_1100 *exp(I*omega1*theta)/(I*omega1)...
    -conj(g1_1100*exp(I*omega1*theta))/(I*omega1)... 
    +g2_1100*exp(I*omega2*theta)/(I*omega2)...
    -conj(g2_1100*exp(I*omega2*theta))/(I*omega2) ...
    +E1100;

w1010 = -g1_1010 *exp(I*omega1*theta)/(I*omega2)...
    -conj(g1_0101*exp(I*omega1*theta))/(I*(2*omega1+omega2))... 
    -g2_1010*exp(I*omega2*theta)/(I*omega1)...
    -conj(g2_0101*exp(I*omega2*theta))/(I*(omega1+2*omega2)) ...
    +E1010*exp(I*(omega1+omega2)*theta);

w1001 = g1_1001 *exp(I*omega1*theta)/(I*omega2)...
    +conj(g1_0110*exp(I*omega1*theta))/(I*(omega2-2*omega1))... 
    +g2_1001*exp(I*omega2*theta)/(I*(2*omega2-omega1))...
    -conj(g2_0110*exp(I*omega2*theta))/(I*omega1) ...
    +E1001*exp(I*(omega1-omega2)*theta);

w0020 = g1_0020 *exp(I*omega1*theta)/(I*(omega1-2*omega2))...
    -conj(g1_0002*exp(I*omega1*theta))/(I*(omega1+2*omega2))... 
    -g2_0020*exp(I*omega2*theta)/(I*omega2)...
    -conj(g2_0002*exp(I*omega2*theta))/(3*I*omega2) ...
    +E0020*exp(2*I*omega2*theta);

w0011 = g1_0011 *exp(I*omega1*theta)/(I*omega1)...
    -conj(g1_0011*exp(I*omega1*theta))/(I*omega1)... 
    +g2_0011*exp(I*omega2*theta)/(I*omega2)...
    -conj(g2_0011*exp(I*omega2*theta))/(I*omega2) ...
    +E0011;

Fq1w1100 = K1*c*(subs(q_1, theta, 0)*subs(diff(w1100, theta), theta, -a1)...
            + subs(w1100, theta, 0)*subs(diff(q_1, theta), theta, -a1)) ...
        + K2*c*(subs(q_1, theta, 0)*subs(diff(w1100, theta), theta, -a2) ...
            + subs(w1100, theta, 0)*subs(diff(q_1, theta), theta, -a2));
        
Fcq1w2000 = K1*c*(subs(conj(q_1), theta, 0)*subs(diff(w2000, theta), theta, -a1)...
            + subs(w2000, theta, 0)*subs(diff(conj(q_1), theta), theta, -a1)) ...
        + K2*c*(subs(conj(q_1), theta, 0)*subs(diff(w2000, theta), theta, -a2) ...
            + subs(w2000, theta, 0)*subs(diff(conj(q_1), theta), theta, -a2));

Fq1w0011 = K1*c*(subs(q_1, theta, 0)*subs(diff(w0011, theta), theta, -a1)...
            + subs(w0011, theta, 0)*subs(diff(q_1, theta), theta, -a1)) ...
        + K2*c*(subs(q_1, theta, 0)*subs(diff(w0011, theta), theta, -a2) ...
            + subs(w0011, theta, 0)*subs(diff(q_1, theta), theta, -a2));

Fq2w1001 = K1*c*(subs(q_2, theta, 0)*subs(diff(w1001, theta), theta, -a1)...
            + subs(w1001, theta, 0)*subs(diff(q_2, theta), theta, -a1)) ...
        + K2*c*(subs(q_2, theta, 0)*subs(diff(w1001, theta), theta, -a2) ...
            + subs(w1001, theta, 0)*subs(diff(q_2, theta), theta, -a2));    
        
Fcq2w1010 = K1*c*(subs(conj(q_2), theta, 0)*subs(diff(w1010, theta), theta, -a1)...
            + subs(w1010, theta, 0)*subs(diff(conj(q_2), theta), theta, -a1)) ...
        + K2*c*(subs(conj(q_2), theta, 0)*subs(diff(w1010, theta), theta, -a2) ...
            + subs(w1010, theta, 0)*subs(diff(conj(q_2), theta), theta, -a2));
        
Fq1cw1001 = K1*c*(subs(q_1, theta, 0)*subs(diff(conj(w1001), theta), theta, -a1)...
            + subs(conj(w1001), theta, 0)*subs(diff(q_1, theta), theta, -a1)) ...
        + K2*c*(subs(q_1, theta, 0)*subs(diff(conj(w1001), theta), theta, -a2) ...
            + subs(conj(w1001), theta, 0)*subs(diff(q_1, theta), theta, -a2));

Fq2w1100 = K1*c*(subs(q_2, theta, 0)*subs(diff(w1100, theta), theta, -a1)...
            + subs(w1100, theta, 0)*subs(diff(q_2, theta), theta, -a1)) ...
        + K2*c*(subs(q_2, theta, 0)*subs(diff(w1100, theta), theta, -a2) ...
            + subs(w1100, theta, 0)*subs(diff(q_2, theta), theta, -a2));      
        
Fcq1w1010 = K1*c*(subs(conj(q_1), theta, 0)*subs(diff(w1010, theta), theta, -a1)...
            + subs(w1010, theta, 0)*subs(diff(conj(q_1), theta), theta, -a1)) ...
        + K2*c*(subs(conj(q_1), theta, 0)*subs(diff(w1010, theta), theta, -a2) ...
            + subs(w1010, theta, 0)*subs(diff(conj(q_1), theta), theta, -a2));

Fq2w0011 = K1*c*(subs(q_2, theta, 0)*subs(diff(w0011, theta), theta, -a1)...
            + subs(w0011, theta, 0)*subs(diff(q_2, theta), theta, -a1)) ...
        + K2*c*(subs(q_2, theta, 0)*subs(diff(w0011, theta), theta, -a2) ...
            + subs(w0011, theta, 0)*subs(diff(q_2, theta), theta, -a2));
    
Fcq2w0020 = K1*c*(subs(conj(q_2), theta, 0)*subs(diff(w0020, theta), theta, -a1)...
            + subs(w0020, theta, 0)*subs(diff(conj(q_2), theta), theta, -a1)) ...
        + K2*c*(subs(conj(q_2), theta, 0)*subs(diff(w0020, theta), theta, -a2) ...
            + subs(w0020, theta, 0)*subs(diff(conj(q_2), theta), theta, -a2));

% The cubic terms, using the formulae in [1].

g1_2100 = conj(D1)*(F1(1,2) + 2*Fq1w1100 + Fcq1w2000);
g2_2100 = conj(D2)*(F1(1,2) + 2*Fq1w1100 + Fcq1w2000);

g1_1011 = conj(D1)*(F1(3,4) + Fq1w0011 + Fq2w1001 + Fcq2w1010);
g2_1011 = conj(D2)*(F1(3,4) + Fq1w0011 + Fq2w1001 + Fcq2w1010);

g1_1110 = conj(D1)*(F1(2,3) + Fq1cw1001 + Fq2w1100 + Fcq1w1010);
g2_1110 = conj(D2)*(F1(2,3) + Fq1cw1001 + Fq2w1100 + Fcq1w1010);

g1_0021 = conj(D1)*(F3(3,4) + 2*Fq2w0011 + Fcq2w0020);
g2_0021 = conj(D2)*(F3(3,4) + 2*Fq2w0011 + Fcq2w0020);

% To compute the normal form we need to evaluate the quadratic and cubic 
% terms at the parameter values for the double hopf bifurcation.

% evaluate g_ijkl
% depending on version of matlab the subs does not/does the convert the 
% symbolic variable to numeric. So an 'eval' is either 
% needed or not. 

% We will need the quadratic terms evaluated in the bifurcation values
try
    
    g1_2000v = eval(subs(g1_2000, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_2000v = eval(subs(g2_2000, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_0020v = eval(subs(g1_0020, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_0020v = eval(subs(g2_0020, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_1100v = eval(subs(g1_1100, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_1100v = eval(subs(g2_1100, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_1010v = eval(subs(g1_1010, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_1010v = eval(subs(g2_1010, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_1001v = eval(subs(g1_1001, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_1001v = eval(subs(g2_1001, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_0200v = eval(subs(g1_0200, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_0200v = eval(subs(g2_0200, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_0002v = eval(subs(g1_0002, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_0002v = eval(subs(g2_0002, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_0011v = eval(subs(g1_0011, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_0011v = eval(subs(g2_0011, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_0101v = eval(subs(g1_0101, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_0101v = eval(subs(g2_0101, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    
    g1_0110v = eval(subs(g1_0110, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
    g2_0110v = eval(subs(g2_0110, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));

catch
    g1_2000v = subs(g1_2000, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_2000v = subs(g2_2000, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_0020v = subs(g1_0020, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_0020v = subs(g2_0020, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_1100v = subs(g1_1100, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_1100v = subs(g2_1100, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_1010v = subs(g1_1010, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_1010v = subs(g2_1010, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_1001v = subs(g1_1001, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_1001v = subs(g2_1001, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_0200v = subs(g1_0200, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_0200v = subs(g2_0200, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_0002v = subs(g1_0002, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_0002v = subs(g2_0002, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_0011v = subs(g1_0011, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_0011v = subs(g2_0011, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_0101v = subs(g1_0101, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_0101v = subs(g2_0101, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    
    g1_0110v = subs(g1_0110, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
    g2_0110v = subs(g2_0110, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
        [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
end

% now the cubic terms
try  
  g1_2100v = eval(subs(g1_2100, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
  g2_2100v = eval(subs(g2_2100, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
          
  g1_1011v = eval(subs(g1_1011, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
  g2_1011v = eval(subs(g2_1011, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));   

  g1_1110v = eval(subs(g1_1110, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
  g2_1110v = eval(subs(g2_1110, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));

  g1_0021v = eval(subs(g1_0021, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
  g2_0021v = eval(subs(g2_0021, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]));
catch
        
  g1_2100v = subs(g1_2100, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
  g2_2100v = subs(g2_2100, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
          
  g1_1011v = subs(g1_1011, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
  g2_1011v = subs(g2_1011, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);   

  g1_1110v = subs(g1_1110, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
  g2_1110v = subs(g2_1110, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);

  g1_0021v = subs(g1_0021, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
  g2_0021v = subs(g2_0021, {K1, K2, omega1, omega2,a1, a2, gamma,c}, ...
              [K1v,K2v, omega1v,omega2v,a1v,a2v, gammav,cv]);
end

%% convert from Wu to Kuznetsov
% quadratic terms
g1_2000v = g1_2000v/2;
g2_2000v = g2_2000v/2;
g1_0020v = g1_0020v/2; 
g2_0020v = g2_0020v/2;
g1_0200v = g1_0200v/2; 
g2_0200v = g2_0200v/2;
g1_0002v = g1_0002v/2; 
g2_0002v = g2_0002v/2;
% cubic terms
g1_2100v = g1_2100v/2;
g2_2100v = g2_2100v/2;
g1_1011v = g1_1011v/1;
g2_1110v = g2_1110v/1;
g1_0021v = g1_0021v/2;
g2_0021v = g2_0021v/2;
% w terms only used to compute manifold do not need conversion
%end of conversions

%% Computation of the cubic terms in the normal form

% To construct the normal form, we first
% transform the four dimensional ODE on the centre manifold after applying
% a near identity transformation to get rid of the quadratic terms. The
% effect of this transformation is that the cubic of the ODE without
% quadratic terms will be polynomial expressions of the quadratic and cubic
% terms. We use the expressions that are previously computed and written in
% [2].

G1_2100v =  g1_2100v + I/omega1v*g1_1100v*g1_2000v ...
    + I/omega2v * (g1_1010v*g2_1100v - g1_1001v*conj(g2_1100v)) ...
    - I/(2*omega1v+omega2v)*(g1_0101v*conj(g2_0200v)) ...
    - I/(2*omega1v-omega2v)*(g1_0110v*g2_2000v);


G1_1011v = g1_1011v + I/omega2v*(g1_1010v*g2_0011v - g1_1001v*conj(g2_0011v)) ...
    +I/omega1v*(2*g1_2000v*g1_0011v - g1_1100v*conj(g1_0011v) ...
    -g1_0011v*conj(g2_0110v) - g1_0011v*g2_1010v) ... % there is a discrepancy here
    -(2*I)/(omega1v+2*omega2v)*g1_0002v*conj(g2_0101v) ...
    -(2*I)/(omega1v-2*omega2v)*(g1_0020v*g2_1001v);

G2_1110v = g2_1110v + I/omega1v*(g1_1100v*g2_1010v - conj(g1_1100v)*g2_0110v)...
    + I/omega2v*(2*g2_0020v*g2_1100v - g2_0011v*conj(g2_1100v) ...
    - g1_1010v*g2_1100v - conj(g1_1001v)*g2_1100v) ...
    + (2*I)/(2*omega1v-omega2v)*g1_0110v*g2_2000v ...
    -(2*I)/(2*omega1v + omega2v)*conj(g1_0101v)*g2_0200v;

G2_0021v = g2_0021v + I/omega2v*g2_0011v*g2_0020v ...
    + I/omega1v*(g1_0011v*g2_1010v - conj(g1_0011v)*g2_0110v) ...
    - I/(2*omega2v+omega1v)*(conj(g1_0002v)*g2_0101v) ...
    - I/(2*omega2v-omega1v)*g1_0020v*g2_1001v;



%% Normal form analysis

% To analyze bifurcations of 2D tori, one has to normalize the fourth- and
% fifth-order terms. The resulting normal form is not unique. If the
% following non-degeneracy conditions hold:
% (HH.0) k*omega1 .neq. l*omega2 for integer k,l>0,k+l .leq. 5 ;
% (HH.1) Re G1_2100(0) .neq. 0 ; These are easy to check
 
disp('The following quatities should all be different from zero.');

condHH0=1e8;
for kit=1:4
for lit=1:5-kit;
    condHH0=min(condHH0,abs(kit*omega1v-lit*omega2v));
end
end
disp(['HH.0: min(abs(k*omega1-l*omega2)) 1\leq k,l, k+l\leq5 is ',num2str(condHH0)])
disp(['HH.1: The real part of G1_2100 is ',num2str(real(G1_2100v))])
% (HH.2) Re g1_1011(0).ne.0 ;
disp(['HH.2: The real part of G1_1011 is ',num2str(real(G1_1011v))])
% (HH.3) Re g2_1110(0) .ne.0 ;
disp(['HH.3: The real part of G2_1110 is ',num2str(real(G2_1110v))])
% (HH.4) Re g2_0021(0) .ne.0 ;
disp(['HH.4: The real part of G2_0021 is ',num2str(real(G2_0021v))])

% (HH.5) the map (kappa_1,kappa_2) maps to (Re lambda_1,Re lambda_2) , 
% where (lambda_1,lambda_2) are eigenvalues of the continuation of the 
% critical equilibrium for small perturbations of (kappa_1,kappa_2)
% from the bifurcation point (kappa_1^*,kappa_2^*) is regular at 
% (kappa_1^*,kappa_2^*) 
%
% To check this note
% lambda + mu + k_1 exp(-lambda a1) + k_2 exp(-lambda a2) = 0
% implies (with all derives partial)
% d lambda/d k_i + exp(-lambda ai) + k_i exp(-lambda ai)[-ai d lambda/d k_i]
%                   + k_j exp(-lambda aj)[-aj d lambda/d k_i] = 0
% implies
%  d lambda/d k_i = -exp(-lambda ai)/(1-ai*k_i exp(-lambda ai)-aj*k_j exp(-lambda aj)
 
dlamdkij=@(ai,aj,ki,kj,lambda) -exp(-lambda*ai)/(1-ai*ki*exp(-lambda*ai)-aj*kj*exp(-lambda*aj)); 

dlamdkap(1,1)=real(dlamdkij(a1v,a2v,K1v,K2v,1i*omega1v));
dlamdkap(1,2)=real(dlamdkij(a2v,a1v,K2v,K1v,1i*omega1v));
dlamdkap(2,1)=real(dlamdkij(a1v,a2v,K1v,K2v,1i*omega2v));
dlamdkap(2,2)=real(dlamdkij(a2v,a1v,K2v,K1v,1i*omega2v));
disp('matrix J used in (kappa1,kapp2)|-->(mu1,mu2) is given by');
disp(dlamdkap);
disp('inv(J)=')
invJ=inv(dlamdkap);
disp(invJ);
disp(' ')
condHH5=det(dlamdkap);
disp(['HH.5: det(dmu/dalpha) is ',num2str(condHH5)])

% The system is now locally orbitally smoothly equivalent near the origin to the complex
% normal form (Gavrilov, 1980)

% Once we have computed all this, the normal form is given by, 
% 
% \dot{r}_1 = r_1(\mu_1 + p_{11}*r_1^2 + p_{12}*r_2^2 + s_1*r_2^4) + O((r_1^2 + r_2^2)^3),
% \dot{r}_2 = r_2(\mu_2 + p_{21}*r_1^2 + p_{22}*r_2^2 + s_2*r_1^4) + O((r_1^2 + r_2^2)^3),
% \dot{\varphi}_1 =  \omega_1 +  O(r_1^2 + r_2^2),
% \dot{\varphi}_2 =  \omega_2 +  O(r_1^2 + r_2^2),
%
% with
%
% {\rm Re}\ P_{11}(0) = {\rm Re}\ G1_2100(0), 
% {\rm Re}\ P_{12}(0) = {\rm Re}\ G1_1011(0)
% {\rm Re}\ P_{21}(0) = {\rm Re}\ G2_1110(0),
% {\rm Re}\ P_{12}(0) = {\rm Re}\ G2_0021(0)
%
% and
%
% p_{kj} = {\rm Re} P_{kj}, s_k={\rm Re} S_k,  k,j=1,2,

p = zeros(2,2);

%p_{11} = Re(P_11} = Re(G_2100)
p(1,1) = real(G1_2100v);

%p_12 = Re(P_12) = Re(G1_1011)
p(1,2) = real(G1_1011v);

%p_21 = Re(P_21) = Re(G2_1110)
p(2,1) = real(G2_1110v);

%p_22 = Re(P_22) = Re(G2_0021)
p(2,2) = real(G2_0021v);



disp(['HH.6: det(p)=',num2str(det(p))])

%% Type of Hopf-Hopf bifrucation

% Quoting [2, section 8.6]
%
% There are essentially different types of bifurcation diagrams of the
% normal form, depending on whether p_11 and p_22 have the same or opposite
% signs. For each of these cases, we also have different subcases.
%
% "Simple" case: p_11 p_22 > 0 
%
% Consider the case p_11<0, p_22 < 0. The case where p_11 and p_22 are
% positive can be reduced to this one by reversion time. Introducing new
% phase space variables and rescaling in time the normal form according to
% xi_1 = -p_11 rho_1, xi_2 = - p_22 rho_2, tau = 2t, yeild
%
% xi_1' = xi_1(\mu_1 - xi_1 - theta xi_2 + Theta xi_2^2)
% xi_2' = xi_2(\mu_2 - delta xi_1 - xi_2 + Delta xi_2^2)
% 
% where theta = p_12/p_22, delta = p_21/p_11, Theta = s_1/p_22^2, 
% Delta = s_2/p_22^2.

% Notice that theta*delta -1 \neq 0 is equivalent to det(p) \neq 0. We can
% assume without loss of generality that theta >= delta. (otherwise reverse
% time and exchange the subscripts in the normal form), the corresponding
% cases follow

%% we first consider case  p(1,1)*p(2,2)>0 considered by Kuznetsov
if p(1,1) > 0 &&  p(2,2) > 0
   disp(['p_11=',num2str(p(1,1)),'>0 & p_22=',num2str(p(2,2)),'>0'])
   disp('Kuznetsov says to reverse time in this case')
   disp('This changes sign of both p & s') 
   disp('Computations that follow are under this transfromation')
   disp('but s has not been computed yet')
   p=-p; % s=-s;
end
   
%% determine bifurcation type
if  p(1,1) < 0 &&  p(2,2) < 0
    disp('Kuznetsov Simple case p11*p22>0')
    %% case considered by Kuznetsov

    theta = p(1,2)/p(2,2);
    delta = p(2,1)/p(1,1);

    if theta - delta < 0
       disp(['theta-delta=',num2str(theta-delta),'<0 Need to reverse time and exchange subscripts of the normal form'])
       theta = p(2,1)/p(1,1);
       delta = p(1,2)/p(2,2);
       disp('Bifurcations below are for swapped parameters')
    elseif theta==delta
       error('theta=delta, not possible, stopping')
    else
       disp(['theta-delta=',num2str(theta-delta),'>0 as required']) 
    end
    
    disp(['theta=',num2str(theta),' delta=',num2str(delta)])
    % Cases of type of H-H bifurcation
    disp('Type of Hopf-Hopf bifurcation:');
    if (theta >0 && delta > 0)
        if (theta*delta >1)
            disp('Case I');
        elseif (theta*delta < 1)
            disp('Case II');
        else
            disp('theta*delta=1; case undetermined');
        end
    elseif (theta > 0 && delta <0)
            disp('Case III');
    elseif (theta <0 && delta < 0)
            if (theta*delta < 1)
                disp('Case IV');
            elseif (theta*delta > 1)
                disp('Case V')
            else
                disp('theta*delta=1; case undetermined');
            end
    else
           error('theta*delta=0, not possible, stopping')
    end
%% we next consider case  p(1,1)*p(2,2)<0 hard case from [2].
elseif p(1,1)*p(2,2)<0
    disp('Kuznetsov Difficult case p11*p22<0')
    disp(['HH.7: p11-p12 is ',num2str(p(1,1)-p(1,2))])
    disp(['HH.8: p21-p22 is ',num2str(p(2,1)-p(2,2))])
    if p(1,1)<0
      disp(['p_11=',num2str(p(1,1)),'>0 & p_22=',num2str(p(2,2)),'>0'])
      disp('reverse time in this case')
      disp('This changes sign of both p & s') 
      disp('Computations that follow are under this transfromation')
      disp('but s has not been computed yet')
      p=-p; % s=-s;
    end
    theta = p(1,2)/p(2,2);
    delta = p(2,1)/p(1,1);    
    disp('Need to compute l_1 in this case: this is NOT done')
    disp('but if theta \geq delta and l_1<0 bifurcations would be:')
    disp(['theta=',num2str(theta),' delta=',num2str(delta)])
    
    if theta - delta < 0
       disp(['theta-delta=',num2str(theta-delta),'<0 Need to reverse time and exchange subscripts of the normal form'])
       theta = p(2,1)/p(1,1);
       delta = p(1,2)/p(2,2);
       disp('Bifurcations below are for swapped parameters')
       disp(['theta=',num2str(theta),' delta=',num2str(delta)])
    end
    
    if theta > delta
        if delta>1 
            disp('difficult Case I');
        elseif theta>1 && delta<1 && theta*delta>1
            disp('difficult Case II');
        elseif delta>0 && theta*delta<1
            disp('difficult Case III');
        elseif delta<0 && theta>0
            disp('difficult Case IV');
        elseif theta<0 && theta*delta<1
            disp('difficult Case V');
        elseif theta<0 && theta*delta>1
            disp('difficult Case VI');
        else
            disp('difficult not elaborated');
        end
    else    
        disp('theta<delta case not elaborated')
    end
else
    disp(['p_11=',num2str(p(1,1)),' p_22=',num2str(p(2,2)),' :Case not yet elaborated'])
end

disp(['theta=',num2str(theta,8),' delta=',num2str(delta,8)])
ttttoc=toc;
disp(['Total Runtime ',num2str(ttttoc),' seconds, or,'])
disp(['Excluding time to find HH point ',num2str(ttttoc+ttoc-tttoc),' seconds'])
