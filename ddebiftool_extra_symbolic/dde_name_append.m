%% append letter to cell array of names
function xnap=dde_name_append(xn,app)
xnap=cell(size(xn));
for i=1:numel(xn)
    xnap{i}=[xn{i},app];
end
end
