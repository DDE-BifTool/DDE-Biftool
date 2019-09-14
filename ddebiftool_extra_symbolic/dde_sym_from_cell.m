%% Convert possibly cells to sym arrays safely
% This wrapper converts if necessary, but does nothing otherwise.
%
% $Id: dde_sym_from_cell.m 309 2018-10-28 19:02:42Z jansieber $
%%
function xa=dde_sym_from_cell(xc)
if iscell(xc)
    xa=cell2sym(xc);
elseif isa(xc,'sym')
    xa=xc;
else
    error('dde_sym_from_cell:input',[...
        'dde_sym_from_cell: input is neither cell (of strings) nor sym array',...
        'but of type %s'],class(xc));
end
