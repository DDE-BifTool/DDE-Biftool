function trimpoint=nmfm_trim_point(point,pref)
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
