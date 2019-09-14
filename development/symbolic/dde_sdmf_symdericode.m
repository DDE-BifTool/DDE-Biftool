function [funcstr,df,v]=dde_sdmf_symdericode(fs,tau,xxs,ps,funcname,varargin)
%% combined functional in equilibrium: derivatives from symbolic expression and matlabcode
%
% code is returned as string containing function
%
% $Id$
%%
%% differentiate expression fs wrt xxs and ps and prepend function itself
fs=fs(:);
tau=tau(:);
[df,v]=dde_sdmf_symdiff(fs,tau,xxs,varargin{:});
df=df([1,1:end]);
%% code generation
funcstr=dde_symcode(df,xxs(:,1),ps,v,sym([]),'funcname',funcname,varargin{:});
end

