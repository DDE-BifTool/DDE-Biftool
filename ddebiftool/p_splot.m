%% Plot spectrum of point with stability information
%%
function p_splot(point,varargin)
% function p_splot(point)
% INPUT:
%	point point whose stability needs plotting
%  'plotaxis' (optional, default gca): axis on which to plot
%
% $Id: p_splot.m 362 2019-07-14 15:49:40Z jansieber $
%
%%
default={'plotaxis',gca};
options=dde_set_options(default,varargin,'pass_on');
if isfield(point,'l0')
    point=struct('kind','stst','stability',point);
elseif isfield(point,'mu')
    point=struct('kind','psol','stability',point);
end    
switch point.kind
    case {'stst','fold','hopf'}
        if isfield(point.stability,'l0')
            root_plt(point.stability.l0,point.stability.l1,[],'plotaxis',options.plotaxis);
        end
        xlabel('\Re(\lambda)');
        ylabel('\Im(\lambda)');
    case 'psol'
        if isfield(point.stability,'mu')
            mult_plt(point.stability.mu,'plotaxis',options.plotaxis);
        end
        xlabel('\Re(\mu)');
        ylabel('\Im(\mu)');
    otherwise
        err=point.kind;
        error('P_SPLOT: point kind %s not recognized.',err);
end
end
