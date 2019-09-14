function funcs=set_rotsymfuncs(func,varargin)
defaults={'rotation',[],'exp_rotation',[]};
[options,pass_on]=dde_set_options(defaults,varargin,'pass_on');
sfuncs=set_symfuncs(func,pass_on{:},'p_vectorized',false);
funcs=set_rotfuncs('sys_rhs',sfuncs.sys_rhs,'sys_dirderi',sfuncs.sys_dirderi,...
    'rotation',options.rotation,'exp_rotation',options.exp_rotation,...
    'x_vectorized',true,'p_vectorized',false,pass_on{:},'sys_cond_reference',true);
end
