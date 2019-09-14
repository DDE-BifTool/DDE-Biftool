function [funcstr,df,v,q]=dde_symdericode(fs,xxs,ps,funcname,varargin)
%% create partial derivatives from symbolic expression and matlabcode
%
% code is returned as string containing function
%
% $Id: dde_symdericode.m 170 2017-03-05 03:39:50Z jansieber $
%%
%% differentiate expression fs wrt xxs and ps and prepend function itself
fs=fs(:);
[df,v,q]=dde_symdiff(fs,xxs,ps,varargin{:});
df=[{fs},df];
%% code generation
funcstr=dde_symcode(df,xxs,ps,v,q,'funcname',funcname,varargin{:});
end