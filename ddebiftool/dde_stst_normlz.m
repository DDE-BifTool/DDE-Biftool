function stst=dde_stst_normlz(stst)
%% normalize steady-state bifurcation points
% (descended from |'stst'| kind)
%
% $Id: dde_stst_normlz.m 366 2019-07-14 21:52:54Z jansieber $
%%
if isfield(stst,'v')
    stst.v=stst.v/norm(stst.v);
end
end
