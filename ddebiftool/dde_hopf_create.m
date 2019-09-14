function [pt,varfields,kind]=dde_hopf_create(varargin)
%% create Hopf point with empty stability and normal form information (unless specified explicitly)
%
% $Id: dde_hopf_create.m 315 2019-01-29 19:42:21Z jansieber $
%%
default={'kind','hopf','parameter',[],'x',[],'v',[],'omega',[],'stability',[],'nmfm',[],'nvec',[],'flag',''};
[pt,dum,userdefined]=dde_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='hopf';
pt=dde_point_overwritefields(pt,userdefined,'stability',[],'nmfm',[],'nvec',[],'flag','');
varfields=struct('x',1,'v',2,'omega',1,'parameter',1);
kind='stst';
end
