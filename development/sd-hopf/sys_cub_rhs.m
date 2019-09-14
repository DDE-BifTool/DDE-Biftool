function f = sys_cub_rhs(xx,par)
%% constant 9 delay DDE approx of epsilon epsilon*u'(t)=-gamma*u(t)-kappa1*u(t-a1-c_1*u(t))-kappa2*u(t-a2-c_2*u(t))
% parameters: par(1)=epsilon, par(2)=gamma, par(3:4)=kappa1:2, par(5:6)=a1:2, par(7:8)=c1:2
epsilon=par(1); gamma=par(2); kappa=par(3:4); c=par(7:8); 

% a2=par(9:11);
% a2=(2a1 a1+a2 2a2);
% a3=par(12:15);
% a3=(3a1 2a1+a2 a1+2a2 3a2);

u=xx(1); ua=xx(2:3);  uaa=xx(4:6); uaaa=xx(7:10);
% u = u(t)
% ua= [u(t-a1) u(t-a2)]
% uaa= [u(t-2a1) u(t-a1-a2) u(t-2a2)]
% uaaa= [u(t-3a1) u(t-2a1-a2) u(t-a1-2a2) u(t-3a2)]

% linear part
f=-gamma*u;
for i=1:2
   % linear part continued
   f=f-kappa(i)*ua(i);
   % quadratic part 
   f=f-gamma*c(i)*u*kappa(i)*ua(i);
   % 33 part
   f=f-0.5*kappa(i)*c(i)^2*u^2*gamma^2*ua(i);
   for j=1:2
      % quadratic part continued
      f=f-c(i)*u*kappa(i)*kappa(j)*uaa(i+j-1);
      % 23 part
      f=f-gamma*c(i)*c(j)*kappa(i)*kappa(j)*u*ua(i)*uaa(i+j-1);
      % 33 part continued
      f=f-kappa(i)*c(i)^2*u^2*gamma*kappa(j)*uaa(i+j-1);
      for m=1:2
        % 23 part continued  
        f=f-c(i)*c(j)*kappa(i)*kappa(j)*u*ua(i)*kappa(m)*uaaa(i+j+m-2);
        % 33part continued
        f=f-0.5*kappa(i)*c(i)^2*u^2*kappa(j)*kappa(m)*uaaa(i+j+m-2);
      end
   end
end
f=f/epsilon;

end
