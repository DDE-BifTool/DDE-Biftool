function [pt,varfields,kind]=dde_psol_create(varargin)
%% create psol point with empty stability information
%
% $Id: dde_psol_create.m 315 2019-01-29 19:42:21Z jansieber $
%%
default={'kind','psol','parameter',[],'mesh',[],'degree',[],'profile',[],...
    'period',[],'stability',[],'nmfm',[],'nvec',[],'flag',''};
[pt,dum,userdefined]=dde_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='psol';
if ~isempty(pt.profile)
    pt=dde_coll_check(pt);
end
pt=dde_point_overwritefields(pt,userdefined,...
    'stability',[],'nmfm',[],'nvec',[],'flag','');
varfields=struct('profile',1,'period',1,'parameter',1);
kind='coll';
end