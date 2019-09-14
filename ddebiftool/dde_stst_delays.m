function taus=dde_stst_delays(funcs,pt,varargin)
%% compute all or requested delays in stst, including first delay=0
%
% $Id: dde_stst_delays.m 299 2018-09-26 09:09:15Z jansieber $
%%
default={'d_nr',[]};
options=dde_set_options(default,varargin,'pass_on');
if ~funcs.tp_del
    taus = pt.parameter(funcs.sys_tau());
    taus = [0, taus]; % First delay zero
    if ~isempty(options.d_nr)
        taus=taus(options.d_nr+1);
    end
else
    %% state-dependent delays
    rg=1:funcs.sys_ntau()+1;
    if ~isempty(options.d_nr)
        rg=options.d_nr+1;
    end
    taus = zeros(1,length(rg)); % First delay zero
    for i = 1:length(rg)
        if rg(i)==1
            continue
        end
        taus(i) = funcs.sys_tau(rg(i)-1,pt.x(:,ones(1,rg(i)-1)),pt.parameter);
    end
end
taus=taus(:);
end
