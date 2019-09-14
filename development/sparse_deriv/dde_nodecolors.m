function coloring=dde_nodecolors(nbhcols)
%% color in grid such that nodes that influence common eqs have different color
% uses greedy graph coloring
%
% input: n x 2 array where each row is a pair of columns that have a
% non-zero entry in a common row
%
% output m x 2 array clr, where clr(:,2) is the list of present columns and
% clr(:,1) is the color of each column
%
%%
col_bd=[1;find(diff(nbhcols(:,1)))+1];
col_bd=[col_bd,[col_bd(2:end)-1;size(nbhcols,1)]];
%% order columns (=nodes) in decreasing degree
[~,ix]=sort(diff(col_bd,[],2),'descend');
ir(ix)=1:length(ix);
nbhreo=reshape(ir(nbhcols(:)),[],2);
nbhreo=sortrows(nbhreo,[1,2]);
reo_bd=[1;find(diff(nbhreo(:,1)))+1];
reo_bd=[reo_bd,[reo_bd(2:end)-1;size(nbhreo,1)]];
columns=nbhreo(reo_bd(:,1),1);
ncolumns=size(columns,1);
clr=NaN(ncolumns,1);
clr(columns(1))=1;
colors=1;
for i=2:ncolumns
    i_nbh=nbhreo(reo_bd(i,1):reo_bd(i,2),2);
    i_nbh=i_nbh(i_nbh<columns(i));
    if isempty(i_nbh)
        clr(columns(i))=1;
        continue;
    else
        used_clr=clr(i_nbh);
        newclr=true(1,length(colors));
        newclr(used_clr)=false;
        if any(newclr)
            clr(columns(i))=colors(find(newclr,1,'first'));
        else
            clr(columns(i))=length(colors)+1;
            colors=[colors,length(colors)+1]; %#ok<AGROW>
        end
    end
end
coloring=unique([clr(columns),ix(nbhreo(reo_bd(:,1),1))],'rows');
end
