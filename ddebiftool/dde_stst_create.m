function [pt,varfields,kind]=dde_stst_create(varargin)
%% create stst point with empty stability and normal form information
%
% $Id: dde_stst_create.m 315 2019-01-29 19:42:21Z jansieber $
%%
default={'kind','stst','parameter',[],'x',[],...
    'stability',[],'nmfm',[],'nvec',[],'flag',''};
[pt,dum,userdefined]=dde_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='stst';
pt=dde_point_overwritefields(pt,userdefined,...
    'stability',[],'nmfm',[],'nvec',[],'flag','');
varfields=struct('x',1,'parameter',1);
kind='';
end