function out=dde_method_check(method,fieldname,mval)
%% check if method has field mfield and if yes compare with mval (if present)
%
% $Id: dde_method_check.m 308 2018-10-28 15:08:12Z jansieber $
%%
out=isfield(method,fieldname);
if ~out
    return
end
mf=method.(fieldname);
if nargin==2
    out=mf;
    return
end
if ischar(mval)
    out=strcmp(mval,mf);
elseif isnumeric(mval)
    out=numel(mf)==numel(mval) && all(mf(:)==mval(:));
end
end
    