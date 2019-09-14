function [hopf_br,augmented] = nmfm_hopf_from_zeho_init(funcs,zeho,radius,freepars,varargin)
%% Initialize branch for continuing the Hopf curve
% going through a zero-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
%
% (adapted by JS)
%
% $Id: nmfm_hopf_from_zeho_init.m 314 2019-01-24 14:28:23Z mmbosschaert $
%%
if ~zeho.nmfm.transcritical 
    %% generic case
    hopf_br=nmfm_eqbif_from_zeho_init(funcs,zeho,radius,freepars,...
        'bifurcation','hopf','evadjust',false);
else
    %% transcritical case: two hopf bifurcations
    g110=zeho.nmfm.g110;
    g200=zeho.nmfm.g200;
    pshift=zeho.nmfm.K*[1,1; 0,real(g110)/g200];
    xshift=[zeros(size(zeho.x)), (-1/g200)*zeho.nvec.q0];
    for i=2:-1:1
        hopf_br(i)=nmfm_eqbif_from_zeho_init(funcs,zeho,radius,freepars,...
        'bifurcation','hopf','evadjust',false,...
        'pshift',pshift(:,i),'xshift',xshift(:,i));
    end
end
augmented=true;
end
