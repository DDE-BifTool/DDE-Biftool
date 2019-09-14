function pt=dde_point_overwritefields(pt,userdefined,varargin)
%% for dde_kind_create functions, overwrite selected fields if not explicitly set
%
% $Id: dde_point_overwritefields.m 368 2019-07-15 07:53:59Z jansieber $
%%
assert(mod(length(varargin),2)==0);
for i=1:2:length(varargin)
    name=varargin{i};
    value=varargin{i+1};
    if ~userdefined.(name)
        pt.(name)=value;
    end
end
end