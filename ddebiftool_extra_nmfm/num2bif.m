function biftype = num2bif(index)
%% Convert an index to a bifurcation type
% Return number of supported types on 'count'
%
% $Id: num2bif.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
if nargin>0
    biftype=bif_num(index,'<-');
else
    biftype=bif_num();
end
end
