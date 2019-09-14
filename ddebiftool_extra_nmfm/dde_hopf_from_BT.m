function hopf=dde_hopf_from_BT(point,varargin)
%% convert Bogdanov-Takens point to (singular!) Hopf
%
% $Id: dde_hopf_from_BT.m 308 2018-10-28 15:08:12Z jansieber $
%%
hopf=dde_hopf_create('point',point,'omega',0,'v',point.nvec.q0);
if isfield(point,'stability')
    hopf.stability=point.stability;
end
end