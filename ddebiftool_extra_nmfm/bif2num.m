function index = bif2num(biftype)
%% Convert bifurcation type to its index
% Return number of supported types on 'count'
%
% $Id: bif2num.m 309 2018-10-28 19:02:42Z jansieber $
%%
if nargin>0
    index=bif_num(biftype,'->');
else
    index=bif_num();
end
end
