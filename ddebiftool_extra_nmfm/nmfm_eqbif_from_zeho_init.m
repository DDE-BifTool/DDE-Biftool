function [bif_br,augmented] = nmfm_eqbif_from_zeho_init(funcs,zeho,radius,freepars,varargin)
%% Initialize branch for continuing the equilibrium bifurcation curve
% going through a zero-Hopf point (fold or Hopf).
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
%
% (adapted by JS)
%
% $Id: nmfm_eqbif_from_zeho_init.m 309 2018-10-28 19:02:42Z jansieber $
%%
default={'bifurcation','','evadjust',true,'pshift',zeho.nmfm.K(:,2),...
    'xshift',zeros(size(zeho.x))};
options=dde_set_options(default,varargin,'pass_on');
p_toeqbif=str2func(['p_to',options.bifurcation]);
biftemplate=p_toeqbif(funcs,zeho);
bif=repmat(biftemplate,1,length(radius));
for i=1:length(radius)
    bif(i).parameter(freepars)=...
        biftemplate.parameter(freepars)+radius(i)*options.pshift';
    bif(i).x=bif(i).x+options.xshift*radius(i);
    if options.evadjust
        bif(i)=p_toeqbif(funcs,bif(i));
    end
end
bif_br=df_brnch(funcs,freepars,options.bifurcation);
bif_br.point=bif;
augmented=true;
end
