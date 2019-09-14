function [pt,vars,kind]=dde_psolfold_create(varargin)
%% create psolfold point with empty stability information (unless specified explicitly)
%
% $Id: dde_psolfold_create.m 340 2019-05-09 19:50:29Z jansieber $
%%
default={'kind','psolfold','parameter',[],'mesh',[],'degree',[],...
    'profile',[],'period',[],...
    'profile_var',[],'parameter_var',[],...
    'parameter_delay',[],'free_par',[],'stability',[]};
[pt,dum,userdefined]=dde_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='psolfold';
pt=dde_point_overwritefields(pt,userdefined,'stability',[],...
    'profile_var',reshape(pt.profile_var,size(pt.profile_var,1),size(pt.profile,2)));
varfields=struct('profile',1,'profile_var',1,'period',1,'parameter',1,...
    'parameter_var',1,'parameter_delay',1);
varfunc=@(p,free_par,base)dde_ind_from_psolbif(p,free_par,base,fieldnames(varfields));
vars={varfields,varfunc};
kind='psol';
end
%%
