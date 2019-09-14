function trimpoint=dde_trim_point(point,pref)
%% Remove all fields of structure point that are not present in pref
% but keep values of point for those fields that ar epresent in pref.
%
% $Id: dde_trim_point.m 309 2018-10-28 19:02:42Z jansieber $
%% 
trimpoint=pref;
field=fieldnames(pref);
for i=1:length(field)
    if isfield(point,field{i})
        trimpoint.(field{i})=point.(field{i});
    else
        trimpoint.(field{i})=[];
    end
end
end
