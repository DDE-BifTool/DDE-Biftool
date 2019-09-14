function [Jp,ircall]=dde_jacobian2pattern(irc,varargin)
%% return struct that helps to obtain sparse 2nd derivative
% with finite differences
%
% input: coumn dimension, index of entries of jacobian: rows,columns of 1st
% dev, cols of 2nd dev
%
% output Jpattern.entries:colors, 
%% find columns that have nonzero entries in common rows
% Acon is incidence matrix of column connection graph
default={'ncol',[]};
options=dde_set_options(default,varargin,'pass_on');
if issparse(irc)
    [ir,ic1]=find(irc);
else
    ir=irc(:,1);
    ic1=irc(:,2);
end
%% combine all possibilities for ic2
irc1=sortrows([ir,ic1],[1,2]);
ir=irc1(:,1);
ic=irc1(:,2);
col_bd=[1;find(diff(ir(:,1)))+1];
col_bd=[col_bd,[col_bd(2:end)-1;size(ir,1)]];
n2c=diff(col_bd,[],2)+1;
n2=n2c.*(n2c+1)/2;
sn2=sum(n2);
ircall=NaN(sn2,3);
k=0;
for i=1:size(col_bd,1)
    [t1,t2]=find(sparse(triu(ones(n2c(i)))));
    icols=ic(col_bd(i,1):col_bd(i,2));
    icols=reshape(icols([t1;t2]),[],2);
    ircall(k+(1:size(icols,1)),1)=ir(col_bd(i,1));
    ircall(k+(1:size(icols,1)),2:3)=icols;
    k=k+size(icols,1);
end
ircall=sortrows(ircall,[1,2,3]);
[irc1,~,lab1]=unique(ircall(:,[1,2]),'rows');
[~,lab2]=ismember(ircall(:,[1,3]),ircall(:,[1,2]),'rows');
Jp.i_sorted=ircall;
Jp1=dde_jacobian1pattern(irc1(:,1:2),'ncol',options.ncol);
Jp.devnumber=Jp1.devnumber(lab1);
Jp.devnumber(:,2)=Jp.devnumber(lab2,1);
%% construct deviations
nlindev=max(Jp.devnumber(:));
[idev1,idev2]=find(sparse(ones(nlindev)));
lindev=(min(1,sparse(Jp.i_sorted(:,2),Jp.devnumber(:,1),...
    ones(size(Jp.devnumber,1),1))));
sel=idev1<=idev2;
dev1=lindev(:,idev1(sel));
dev2=lindev(:,idev2(sel));
Jp.dev=[dev1+dev2,dev1-dev2]/2;
[~,lmap]=ismember([idev2(~sel),idev1(~sel)],[idev1,idev2],'rows');
Jp.expand=@(fvals)jac_expand(options.ncol,fvals,Jp,lmap,sel,idev1,idev2);
end
%%
function df=jac_expand(nx,fdiffs,pat,lmap,sel,idev1,idev2)
[nf,ndev]=size(fdiffs);
fdiffs=fdiffs(:,1:ndev/2)-fdiffs(:,ndev/2+1:end);
fdiffall(:,sel)=fdiffs;
fdiffall(:,~sel)=fdiffall(:,lmap);
[~,lin_indev]=ismember(pat.devnumber,[idev1,idev2],'rows');
%% insert values into sparse Jacobian
df.ind=pat.i_sorted;
lin_ind=sub2ind(size(fdiffall),pat.i_sorted(:,1),lin_indev);
df.vals=fdiffall(lin_ind);
%% remove exact zero entries
isnz=df.vals~=0;
df.vals=df.vals(isnz);
df.ind=df.ind(isnz,:);
%% duplicate symmetric entries
idiff=df.ind(:,2)<df.ind(:,3);
df.ind=[df.ind;df.ind(idiff,[1,3,2])];
df.vals=[df.vals;df.vals(idiff)];
[df.ind,ix]=sortrows(df.ind);
df.vals=df.vals(ix);
df.size=[nf,nx,nx];
end
