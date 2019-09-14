function [funcstr,df,v,q]=dde_symdericode(fs,xxs,ps,funcname,varargin)
%% create partial derivatives from symbolic expression and matlabcode
%
% code is returned as string containing function
%
% $Id: dde_symdericode.m 309 2018-10-28 19:02:42Z jansieber $
%%
%% differentiate expression fs wrt xxs and ps and prepend function itself
fs=fs(:);
[df,v,q]=dde_symdiff(fs,xxs,ps,varargin{:});
df=[{fs},df];
%% code generation
funcstr=dde_symcode(df,xxs,ps,v,q,'funcname',funcname,varargin{:});
end