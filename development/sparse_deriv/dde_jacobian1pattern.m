function Jpattern=dde_jacobian1pattern(irc,varargin)
%% return struct that helps to obtain sparse Jacobian with finite differences
%
%% Input: 
% 
% * irc: sparse matrix of which the index set ([ir,ic]=find(in1) is the
% fill pattern, or row and columns indices of fill pattern (in mx2 form)
%
%% Output Jpattern fields
%
% * i_sorted: m x 2 indices (rows, columns)
% * devnumber: number of deviation that will recover entry at i_sorted 
% * dev: deviations (sparse, vectors of 0 and 1's)
% * expand: function. After evaluating Jdev(:,i)=df(x)dev(:,i) for all i,
% call J=Jpattern.expand(Jdev) to obtain full Jacobian.
%
% The function dde_nodecolors partitions the column index set 1:nx into
% groups such that columns from each group have no common row indices in
% the fill pattern. This is done by a greedy graph coloring after sorting
% the nodes for decreasing degree. This algorithm is known to give
% solutions that have no worse than twice the minimal number of groups.
%% Process inputs
% 2nd/3rd input is incidence matrix of column connection graph
default={'dev',[],'ncol',[]};
options=dde_set_options(default,varargin,'pass_on');
if issparse(irc)
    [ir,ic]=find(irc);
else
    ir=irc(:,1);
    ic=irc(:,2);
end
if isempty(options.ncol)
    if issparse(irc)
        nx=size(irc,2);
    elseif ~isempty(options.dev)
        nx=size(options.dev,1);
    else
        nx=max(ic);
    end
else
    nx=options.ncol;
end
%% label column/rows skipping empty rows, columns
[colnames,~,ic]=unique(ic);
[rownames,~,ir]=unique(ir);
%% find out which columns have common row entry patterns
Ainc=sparse(ir,ic,ones(size(ir)));
[nbhcol1,nbhcol2]=find(Ainc'*Ainc);
nbhcols=sortrows([nbhcol1,nbhcol2]); % pairs of "connected" cols
%% partition connection graph into clusters with non-neighboring columns
colors2columns=dde_nodecolors(nbhcols);
colors=colors2columns(:,1);
columns=colors2columns(:,2);
%% sort entries according to colored columns
ncols=length(columns);
ix(columns)=1:ncols;
[~,kx]=sort(ix(ic));
ir=ir(kx);
ic=ic(kx);
column_bd=[1;find(diff(ic))+1];
column_bd=[column_bd,[column_bd(2:end)-1;length(ic)]];
irc=[rownames(ir),colnames(ic)];
[ircsort,ix]=sortrows(irc);
clr=NaN(length(ic),1);
for i=1:size(column_bd,1)
    clr(column_bd(i,1):column_bd(i,2))=colors(i);
end
clr=clr(ix);
%% dev are the deviations needed for the Jacobian
dev1=(min(1,sparse(ircsort(:,2),clr,ones(length(clr),1),nx,max(clr))));
if ~isempty(options.dev)
    dev=[dev1+options.dev,dev1-options.dev]*0.5;
    order=2;
else
    order=1;
    dev=dev1;
end
Jpattern=struct('i_sorted',ircsort,'devnumber',clr,...
    'dev',dev,'order',order);
Jpattern.expand=@(fvals)jac_expand(nx,fvals,Jpattern);
end
%%
function df=jac_expand(nx,fdiffs,pat)
[nf,ndev]=size(fdiffs);
if pat.order==2
    fdiffs=fdiffs(:,1:ndev/2)-fdiffs(:,ndev/2+1:end);
    ndev=ndev/2;
end
df=sparse(nf,nx);
ldf_ind=sub2ind(size(df),pat.i_sorted(:,1),pat.i_sorted(:,2));
ldev_ind=sub2ind([nf,ndev],pat.i_sorted(:,1),pat.devnumber);
df(ldf_ind)=fdiffs(ldev_ind);
end
%%
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
