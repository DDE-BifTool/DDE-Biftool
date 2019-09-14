function tau=poscontrol_tau(delay_nr,xx,par)

if delay_nr==1
  tau=par(1)+0*xx(1,1,:);
elseif delay_nr==2
  tau=xx(2,1,:);
elseif delay_nr==3
  tau=par(1)+xx(2,1,:);
end;
end
