function [stst2,stpcond]=dde_stst_from_stst(stst,varargin)
%% branch off at stst which is approximate branch point (not corrected)
%
% optional arguments: 
%
% * 'funcs': r.h.s. functions structure to compute char. matrix if needed
%
% $Id: dde_stst_from_stst.m 308 2018-10-28 15:08:12Z jansieber $
%%
default={'radius',1e-3};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
stst_ext=dde_fold_from_stst(stst,pass_on{:});
stst2=stst;
stst2.x=stst2.x+options.radius*stst_ext.v;
stpcond=dde_stst_create('x',stst_ext.v,'parameter',stst.parameter*0);
end
