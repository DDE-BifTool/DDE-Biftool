function [fold_br,augmented] = nmfm_fold_from_zeho_init(funcs,zeho,radius,freepars,varargin)
%% Initialize branch for continuing the fold curve
% going through a zero-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
%
% (adapted by JS)
%
% $Id: nmfm_fold_from_zeho_init.m 314 2019-01-24 14:28:23Z mmbosschaert $
%%
if ~zeho.nmfm.transcritical 
    %% generic case
[fold_br,augmented]=nmfm_eqbif_from_zeho_init(funcs,zeho,radius,freepars,...
    'bifurcation','fold','evadjust',true);
else
    %% transcritical case: two hopf bifurcations
    pshift=zeho.nmfm.K(:,2);
    [fold_br,augmented]=nmfm_eqbif_from_zeho_init(funcs,zeho,radius,freepars,...
    'bifurcation','fold','evadjust',true,'pshift',pshift);
    fold_br.tangent=@(p) [];
end
end
