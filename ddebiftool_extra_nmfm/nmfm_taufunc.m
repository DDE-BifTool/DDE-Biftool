function tau=nmfm_taufunc(order,it,xx,p,dp,itau) %#ok<INUSL>
%% when converting DDE with parameter delays to sd-DDE, 
% this is the sys_tau and sys_dirdtau
%
% $Id: nmfm_taufunc.m 309 2018-10-28 19:02:42Z jansieber $
%%
if order==0
    tau=p(1,itau(it),:);
elseif order==1
    tau=dp(1,itau(it),:);
else
    tau=0*dp(1,1,:);
end
end
