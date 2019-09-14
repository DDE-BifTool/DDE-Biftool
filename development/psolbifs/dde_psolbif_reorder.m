function ind=dde_psolbif_reorder(ind,fnames)
ic=cellfun(@(name)ind.(name),fnames,'uniformoutput',false);
iar=reshape(cat(2,ic{:}),size(ic{1},1),length(fnames),[]);
iarswap=permute(iar,[1,3,2]);
for i=1:length(fnames)
    ind.(fnames{i})=iarswap(:,:,i);
end
end
