function grid=dde_coll_set_grid(purpose,degree,varargin)
%% define various node grids for polynomial interpolation
% $Id: dde_coll_set_grid.m 369 2019-08-27 00:07:02Z jansieber $
%%
default={'type','equidistant','subset',@(x)x(2:end)};
options=dde_set_options(default,varargin,'pass_on');
if isa(options.type,'function_handle')
    grid=options.type(purpose,degree);
    return
elseif isnumeric(options.type) && ~isempty(options.type)
    grid=options.type;
    return
end
if isempty(options.type) 
    if ~strcmp(purpose,'collocation')
        options.type='equidistant';
    else
        options.type='gauss';
    end
end
switch options.type
    case {'gauss','legendre'}
        grid=legpts(degree,[0,1]);
        grid=grid(:).';
    case {'cheb','chebyshev'}
        grid=dde_coll_chebxwt(degree,1:degree+1,[0,1]);
    case {'linear','equidistant'}
        grid=linspace(0,1,degree+1);
end
switch purpose
    case 'collocation'
        if length(grid)>degree
            grid=options.subset(grid);
        end
end
end
