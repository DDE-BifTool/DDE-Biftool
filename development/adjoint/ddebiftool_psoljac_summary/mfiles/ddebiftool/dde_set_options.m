function [options,passed_on]=dde_set_options(defaults,userargs,pass_on)
%% parses varargin and assigns fields of structure options
% unknown arguments are passed on into cell array passed_on if pass_on is
% present and non-empty or false, otherwise and error message is generated
%
% $Id: dde_set_options.m 142 2017-02-08 12:16:41Z jansieber $
%
passed_on={};
% wrap cell arguments to avoid generating multiple structs
if isstruct(defaults)
    options=defaults;
elseif iscell(defaults)
    for i=1:length(defaults)
        if iscell(defaults{i})
            defaults{i}=defaults(i);
        end
    end
    options=struct(defaults{:});
else
    error('defaults not recognized\n');
end
if nargin<3 || isempty(pass_on)
    pass_on=false;
end
if length(userargs)~=1
    for i=1:2:length(userargs)
        if isfield(options,userargs{i})
            options.(userargs{i})=userargs{i+1};
        else
            if ~pass_on
                error('option ''%s'' not recognized\n',userargs{i});
            else
                passed_on={passed_on{:},userargs{i},userargs{i+1}}; %#ok<CCAT>
            end
        end
    end
else
    userargs=userargs{1};
    if ~isstruct(userargs)
        error('option ''%s'' not recognized\n',userargs{i});
    end
    passed_on={};
    fnames=fieldnames(userargs);
    for i=1:length(fnames)
        if isfield(options,fnames{i})
            options.(fnames{i})=userargs.(fnames{i});
        else
            passed_on={passed_on{:},fnames{i},userargs.(fnames{i})}; %#ok<CCAT>
        end
    end
end
end
