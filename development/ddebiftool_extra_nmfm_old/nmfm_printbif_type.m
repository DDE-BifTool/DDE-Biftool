function nmfm_printbif_type(str,occurence)
%% shortcut printing out number of occurences of type str
%
% $Id: nmfm_printbif_type.m 109 2015-08-31 23:45:11Z jansieber $
%
if sum(occurence)>0
    fprintf('%d %s ',sum(occurence),str);
end
end
