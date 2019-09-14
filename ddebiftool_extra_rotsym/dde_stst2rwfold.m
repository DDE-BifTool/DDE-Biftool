function data=dde_stst2rwfold(funcs,info,varargin)
[pffuncs,info]=SetupRWFold(funcs,info,1,'correc',false,'dir',[],varargin{:});
data.info=info;
data.funcs=pffuncs;
end