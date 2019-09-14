function taus=dde_coll_delays(funcs,pt,varargin)
%% compute all or requested delays in stst, including first delay=0
%
% $Id: dde_coll_delays.m 362 2019-07-14 15:49:40Z jansieber $
%%
default={'d_nr',[],'t',pt.mesh};
options=dde_set_options(default,varargin,'pass_on');
coll=dde_coll_map(funcs,pt,'c_is_tvals',true,'c',options.t,'nderivs',0);
taus=coll.tau';
if ~isempty(options.d_nr)
    taus=taus(options.d_nr+1,:);
end
end
