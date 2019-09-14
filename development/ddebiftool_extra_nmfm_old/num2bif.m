function biftype = num2bif(index)
%% Convert an index to a bifurcation type
% Return number of supported types on 'count'
%
% $Id: num2bif.m 112 2015-09-02 00:31:41Z jansieber $
%
%%
if nargin>0
    biftype=bif_num(index,'<-');
else
    biftype=bif_num();
end
end
