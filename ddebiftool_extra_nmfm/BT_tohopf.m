function hopf=BT_tohopf(funcs,btpoint,freqs) %#ok<INUSD,INUSL>
%% convert Bogdanov-takens point into a Hopf point (wrapper)
%
% $Id: BT_tohopf.m 309 2018-10-28 19:02:42Z jansieber $
%%
hopf=dde_hopf_from_BT(btpoint);
end 
  