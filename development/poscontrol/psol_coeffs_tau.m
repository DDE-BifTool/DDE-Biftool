function [y,W,neqs]=psol_coeffs_tau(tctau,pt,varargin)
%% Returns matrix W, W' (and W'') mapping base points to collocation points 
% W{i}{j} is mapping of base points to (i-1)th derivative at time t-tau(j)
% where counting of tau starts with delay tau(1)=0.
%% optional parameters
% c pertmis t ospecify collocation points inside collocation interval,
% scaled to [0,1] (if c_is_tvals is false). If c_is_tvals is true then the
% vector c is treated as the requested time points on the entire interval.
% 'kron' switches on/off if W should be expanded to system size, 'nderivs'
% sets whether W, W' and W'' should be computed.
default={'c',[],'c_is_tvals',false,'nderivs',1,'wrapJ',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
for j=options.nderivs:-1:0
    [y{j+1},W{j+1}]=psol_eva_jac(pt.profile,pt.mesh,c_tau,pt.degree,pass_on{:},'wrapJ',options.wrapJ,'diff',j);
end
end

