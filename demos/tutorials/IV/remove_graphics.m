function remove_graphics()
s=evalin('base','whos');
sc={s.class};
id='matlab.graphics';
for i=1:length(sc)
    if strncmp(sc{i},id,length(id))
        evalin('base',['clear ',s(i).name]);
    end
end
end
