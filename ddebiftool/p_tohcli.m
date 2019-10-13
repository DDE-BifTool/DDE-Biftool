function [hcli,stst]=p_tohcli(point,varargin)
%% convert point to connecting orbit
% INPUT:
%     point a periodic solution near a homoclinic solution
%           alternatively an initial point in a hcli structure,
%           where a good starting guess for the profile and steady
%           states are available
%     named (but mandatory) 'funcs': problem functions
% OUTPUT:
%     hcli a starting value to compute the exact homoclinic or
%     heteroclinic solution  
%
% (c) DDE-BIFTOOL v. 2.02, 16/6/2002
%
% $Id: p_tohcli.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
%% for backward compatibility re-organize arguments
args=varargin;
if ~isfield(point,'kind')
    funcs=point;
    point=args{1};
    args=[args(2:end),{'funcs',funcs}];
end
[hcli,stst]=dde_apply({'dde_hcli_from_',point.kind,''},point,args{:});
end
