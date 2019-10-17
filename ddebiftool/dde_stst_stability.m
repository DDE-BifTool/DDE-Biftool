function stab=dde_stst_stability(varargin)
%% create stst stability field
%
% $Id$
%%
default={'h',[],'l0',NaN(1,0),'l1',NaN(0,1),'err',NaN(1,0),...
    'v',NaN(1,0),'w',NaN(1,0),'discarded',[]};
stab=dde_set_options(default,varargin,'pass_on','stability');
end
