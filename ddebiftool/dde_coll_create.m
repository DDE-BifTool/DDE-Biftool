function [pt,varfields,kind]=dde_coll_create(varargin)
%% create psol point with empty stability information
%
% $Id: dde_coll_create.m 348 2019-06-19 13:09:01Z jansieber $
%%
default={'kind','coll','parameter',[],'mesh',[],'degree',[],'profile',[],...
    'period',[],'stability',[]};
pt=dde_set_options(default,varargin,'pass_on','point');
pt.kind='coll';
if ~isempty(pt.profile)
    pt=dde_coll_check(pt);
end
varfields=struct('profile',1,'period',1,'parameter',1);
kind='';
end