function tau_eva=p_tau(funcs,point,d_nr,t)
%% evaluate (state-dependent) delays along orbit
% function [tau_eva]=p_tau(point,d_nr,t)
% INPUT:
%       point a point
%       d_nr number(s) of delay(s) to evaluate
%       t (optional) point(s) where to evaluate
% OUTPUT:
%       tau_eva value(s) of evaluated delay(s)
%
% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: p_tau.m 362 2019-07-14 15:49:40Z jansieber $
%
if nargin<4
    t_args={};
else
    t_args={'t',t};
end
tau_eva=dde_apply({'dde_',point.kind,'_delays'},funcs,point,'d_nr',d_nr,t_args{:});
end
