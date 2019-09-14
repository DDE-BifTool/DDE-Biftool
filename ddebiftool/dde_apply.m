function varargout=dde_apply(fcn,varargin)
%% wrapper around feval to permit easy overloading for extended structs
fname0=[fcn{:}];
fname=fname0;
fun_ex=exist(fname,'file');
kind=fcn{2};
while ~fun_ex
    [~,~,kind]=feval(['dde_',kind,'_create']);
    if isempty(kind)
        error('dde_apply:function','dde_apply: %s cannot be applied',fname0);
    end
    fname=[fcn{1},kind,fcn{3}];
    fun_ex=exist(fname,'file');
end
varargout=cell(1,nargout);
[varargout{:}]=feval(fname,varargin{:});
end
