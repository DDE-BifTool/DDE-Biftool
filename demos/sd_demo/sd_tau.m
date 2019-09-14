%% User-provided state-dependent delays
% Implementation for tutorial sd_demo:
%
%   function dtau=sd_tau(k,xx,par)
%
% returns delay $\tau_k$, depending on |xx(:,1:k)|
% ($[x(t),\ldots,x(t-\tau_{k-1})]$) and |par|.
%
% $Id: sd_tau.m 20 2014-04-11 19:27:33Z jan.sieber $
%
%%
function tau=sd_tau(k,xx,par)
pad=zeros(1,size(xx,3));
if k==1
    tau=par(10)+pad;
elseif k==2
    tau=par(11)+pad;
elseif k==3
    tau=2+par(5)*par(10)*xx(2,1,:).*xx(2,2,:);
elseif k==4
    tau=1-1./(1+xx(2,3,:).*xx(1,1,:));
elseif k==5
    tau=xx(4,1,:);
elseif k==6
    tau=xx(5,1,:);
end
end

