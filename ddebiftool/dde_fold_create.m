function [pt,varfields,kind]=dde_fold_create(varargin)
%% create fold point with empty stability and normal form information
%
% $Id: dde_fold_create.m 315 2019-01-29 19:42:21Z jansieber $
%%
default={'kind','fold','parameter',[],'x',[],'v',[],'stability',[],'nmfm',[],'nvec',[],'flag',''};
[pt,dum,userdefined]=dde_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='fold';
pt=dde_point_overwritefields(pt,userdefined,'stability',[],'nmfm',[],'nvec',[],'flag','');
varfields=struct('x',1,'v',1,'parameter',1);
kind='stst';
end
