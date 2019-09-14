function data=dde_psol2mwPD(funcs,info,varargin)
data=dde_psol2PD(funcs,info,'nremove',[1,1],varargin{:});
end