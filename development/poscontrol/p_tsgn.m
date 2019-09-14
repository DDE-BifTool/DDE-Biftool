function [delay_nr,tz]=p_tsgn(funcs,point)
%% find negative delays
% function [delay_nr,tz]=p_tsgn(point)
% INPUT:
%       funcs problem function
%       point a point 
% OUTPUT:
%	delay_nr number of negative delay 
%       tz (for psol only) we want tau(tz)=0 and dtau/dt(tz)=0
%
% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: p_tsgn.m 19 2014-04-11 14:15:36Z jan.sieber $
%
%%
d=funcs.sys_ntau(); % number of delays 
refine_levels=5;
delay_nr=[];
tz=[];
tau_all=p_tau(funcs,point,1:d);
switch point.kind
    case {'stst','hopf','fold'}
        delay_nr=find(tau_all<0,1,'first');
    case 'psol'
        %% interpolate tau_eva on refined mesh
        mesh=point.mesh;
        for i=1:refine_levels
            mesh=[mesh(1:end-1);(mesh(1:end-1)+mesh(2:end))*0.5];
            mesh=[mesh(:)',1];
        end
        tau_all_ref=psol_eva(setfield(point,'profile',tau_all),mesh); %#ok<SFLD>
        % check sign of delay
        [tau_n, ind_tau]=min(tau_all_ref,[],2);
        tau_neg=find(tau_n<0,1,'first');
        if ~isempty(tau_neg)
            tz=mesh(ind_tau(tau_neg));
            delay_nr=tau_neg;
        end
    case 'hcli'
        error('P_TSGN: this routine is not (yet) implemented for connecting orbits.');
    otherwise
        error('P_TSGN: point kind %s not recognized.',point.kind);
end
end
