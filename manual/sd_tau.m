function tau=sd_tau(ind,xx,par)
if ind==1
  tau=par(10);
elseif ind==2
  tau=par(11);
elseif ind==3
  tau=2+par(5)*par(10)*xx(2,1)*xx(2,2);
elseif ind==4
 tau=1-1/(1+xx(2,3)*xx(1,1));
elseif ind==5
 tau=xx(4,1);
elseif ind==6
 tau=xx(5,1);
end
end
