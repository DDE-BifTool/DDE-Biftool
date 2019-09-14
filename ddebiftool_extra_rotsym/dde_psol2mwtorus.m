function data=dde_psol2mwtorus(funcs,info,varargin)
data=dde_psol2torus(funcs,info,'nremove',[1,1],varargin{:});
end
