function hopf=dde_hopf_from_zeho(zeho,varargin)
%% extract hopf point from zeho
%
% $Id: dde_hopf_from_zeho.m 308 2018-10-28 15:08:12Z jansieber $
%%
if isempty(zeho.nvec)
    hopf=dde_hopf_from_stst(zeho,varargin{:});
else
    hopf=dde_hopf_create('point',zeho,'v',zeho.nvec.q1,'omega',zeho.nvec.omega,varargin{:});
end
end

