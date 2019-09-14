function data=dde_psol2mwPOfold(funcs,info,varargin)
[pffuncs,info]=SetupMWFold(funcs,info,1,'correc',false,'dir',[],varargin{:});
data.info=info;
data.funcs=pffuncs;
end