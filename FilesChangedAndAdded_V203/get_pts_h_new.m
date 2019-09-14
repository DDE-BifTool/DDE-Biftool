function pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,points)
  
% function pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,eigA0)
% INPUT:
%	AA,tau,real_min,delta_real_min,nb_nu,eigA0
% OUTPUT:
%	pts
%   
% COMMENT: is implemented for m=1, 2, 3.

% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007
  
theta=linspace(0,2*pi,nb_nu+1);
theta=theta(1:end-1);
nu=complex(sin(theta),cos(theta));

r_min_d=real_min-delta_real_min;
ss=exp(-real_min*tau);

E=hlp_get_pts_h_new(AA,nu,ss);
points=[points,E];

rp=real(points);
ip=imag(points);
inds_p=find(rp>=r_min_d);
rp=rp(inds_p);
ip=ip(inds_p);
%

pts=complex(rp,ip);
%
%% ALTERNATIVE:
%% Note: convhull's option 'Pp' does not always 
%% work in the case that "the initial hull is narrow" ...
%% rp=[rp,rp];
%% ip=[ip,-ip];
%% %
%% K=convhull(rp,ip,{'Qt','Pp'});
%% conv_hull=complex(rp(K),ip(K));
%% %
%% pts=conv_hull(2:end);
%% % Could throw away half of these points actually ...
%

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E=hlp_get_pts_h_new(CC,nu,ss)

m=length(CC)-1;  
nb_nu=length(nu);

E=[];

% m ...
switch m
 case 1
  for i2=1:ceil((nb_nu+1)/2)
    e=eig(CC{1}+ss(1)*nu(i2)*CC{2});
    E=[E,e.'];
  end
 case 2
  for i2=1:nb_nu
    for i3=1:ceil((nb_nu+1)/2)
      e=eig(CC{1}+ss(1)*nu(i2)*CC{2}+ss(2)*nu(i3)*CC{3});
      E=[E,e.'];
    end
  end        
 case 3
  for i2=1:nb_nu
    for i3=1:nb_nu
      for i4=1:ceil((nb_nu+1)/2)
	e=eig(CC{1}+ss(1)*nu(i2)*CC{2}+ss(2)*nu(i3)*CC{3}+ss(3)*nu(i4)*CC{4});
	E=[E,e.'];
      end
    end
  end     
end

return;
