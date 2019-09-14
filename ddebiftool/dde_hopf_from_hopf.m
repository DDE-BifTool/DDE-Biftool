function hopf=dde_hopf_from_hopf(hopf,varargin)
%% wrapper for dde_hopf_from_stst, switching to different frequency
%
% $Id: dde_hopf_from_hopf.m 308 2018-10-28 15:08:12Z jansieber $
%%
hopf=dde_hopf_from_stst(hopf,varargin{:});
end
