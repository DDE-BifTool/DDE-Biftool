function nmfm_printbif_type(str,occurence)
%% shortcut printing out number of occurences of type str
%
% $Id: nmfm_printbif_type.m 309 2018-10-28 19:02:42Z jansieber $
%
if sum(occurence)>0
    fprintf('%d %s ',sum(occurence),str);
end
end
