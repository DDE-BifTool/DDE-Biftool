function index = bif2num(biftype)
%% Convert bifurcation type to its index
% Return number of supported types on 'count'
%
% $Id: bif2num.m 112 2015-09-02 00:31:41Z jansieber $
%%
if nargin>0
    index=bif_num(biftype,'->');
else
    index=bif_num();
end
end
