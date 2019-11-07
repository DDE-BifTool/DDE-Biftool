function trimpoint=dde_trim_point(point,pref)
%% Empty all fields of structure point that are not present in pref
% but keep values of point for those fields that ar epresent in pref.
%
% $Id$
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
