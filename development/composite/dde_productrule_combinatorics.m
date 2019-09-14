function cf=dde_productrule_combinatorics(maxorder,varargin)
default={'xorder',maxorder,'yorder',maxorder,'sort','x'};
options=dde_set_options(default,varargin,'pass_on');
%% differentiate g(x(s),y(s)) 1..maxorder (1..p) times
% Results is sum of i*D[j,k](g){D[1]x^(m(1),...D[p]x^m(p), D[1]y^(n(1),...D[p]y^n(p)}
% coefficient for qth derivativeis stored at
% i=cf(q).fac
% j=cf(q).gx;
% k=cf(q).gy;
% m(r)=cf(q).xo;
% n(r)=cf(q).yo;
cf=repmat(struct('fac',[1,1],'gx',[1,0],'gy',[0,1],...
    'xo',[[1,0];zeros(maxorder-1,2)],'yo',[[0,1];zeros(maxorder-1,2)]),maxorder,1);
for i=1:maxorder-1
    nc=cf(i);
    numc=length(nc.fac);
    nrnew=sum(nc.xo(:)>0)+sum(nc.yo(:)>0)+2*numc;
    xo=NaN(maxorder,nrnew);
    yo=NaN(maxorder,nrnew);
    fac=NaN(1,nrnew);
    gx=NaN(1,nrnew);
    gy=NaN(1,nrnew);
    ind=0;
    for j=1:numc
        %% differentiate g wrt x, add 1st deriv of x
        ind=ind+1;
        fac(ind)=nc.fac(j);
        gx(ind)=nc.gx(j)+1;
        gy(ind)=nc.gy(j);
        xo(:,ind)=nc.xo(:,j);
        xo(1,ind)=xo(1,ind)+1;
        yo(:,ind)=nc.yo(:,j);
        %% differentiate g wrt y, add 1st deriv of y
        ind=ind+1;
        fac(ind)=nc.fac(j);
        gx(ind)=nc.gx(j);
        gy(ind)=nc.gy(j)+1;
        xo(:,ind)=nc.xo(:,j);
        yo(:,ind)=nc.yo(:,j);
        yo(1,ind)=yo(1,ind)+1;
        ix=find(nc.xo(:,j)>0);
        iy=find(nc.yo(:,j)>0);
        for k=1:length(ix)
            ind=ind+1;
            fac(ind)=nc.fac(j);
            gx(ind)=nc.gx(j);
            gy(ind)=nc.gy(j);
            xo(:,ind)=nc.xo(:,j);
            yo(:,ind)=nc.yo(:,j);
            %% differentiate any derivative of x once more
            xo(ix(k)+1,ind)=xo(ix(k)+1,ind)+1;
            fac(ind)=fac(ind)*xo(ix(k),ind);
            xo(ix(k),ind)=xo(ix(k),ind)-1;
        end
        for k=1:length(iy)
            ind=ind+1;
            fac(ind)=nc.fac(j);
            gx(ind)=nc.gx(j);
            gy(ind)=nc.gy(j);
            xo(:,ind)=nc.xo(:,j);
            yo(:,ind)=nc.yo(:,j);
            %% differentiate any derivative of x once more
            yo(iy(k)+1,ind)=yo(iy(k)+1,ind)+1;
            fac(ind)=fac(ind)*yo(iy(k),ind);
            yo(iy(k),ind)=yo(iy(k),ind)-1;
        end
    end
    sxo=3;
    syo=sxo+maxorder;
    %% check for equal entries
    mat=[gx;gy;xo;yo;fac]';
    mat=sortrows(mat);
    [~,iu]=unique(mat(:,1:end-1),'rows');
    iu(end+1)=size(mat,1)+1; %#ok<AGROW>
    nfac=zeros(1,length(iu)-1);
    for j=1:length(iu)-1
        nfac(j)=sum(mat(iu(j):iu(j+1)-1,end));
    end
    mat=[mat(iu(1:end-1),1:end-1),nfac'];
    switch options.sort
        case 'x'
            mat=sortrows(mat,sxo+(maxorder-1:-1:0));
        case 'y'
            mat=sortrows(mat,syo+(maxorder-1:-1:0));
    end
    mat=mat(end:-1:1,:);
    mat=mat';
    cf(i+1).fac=mat(end,:);
    cf(i+1).gx=mat(1,:);
    cf(i+1).gy=mat(2,:);
    cf(i+1).xo=mat(sxo+(0:maxorder-1),:);
    cf(i+1).yo=mat(syo+(0:maxorder-1),:);
end
cf=cleanup(cf,'x',options.xorder);
cf=cleanup(cf,'y',options.yorder);
for i=1:maxorder
    if isfield(cf(i),'xo')
        cf(i).xo=cf(i).xo(1:i,:);
    end
    if isfield(cf(i),'yo')
        cf(i).yo=cf(i).yo(1:i,:);
    end
end
end
function cf=cleanup(cf,fname,order)
for i=1:length(cf)
    cfi=cf(i);
    keep=~any(cfi.([fname,'o'])(order+1:end,:)>0);% & cfi.(['g',fname])==1;
    cfi.fac=cfi.fac(keep);
    cfi.gx=cfi.gx(keep);
    cfi.gy=cfi.gy(keep);
    cfi.xo=cfi.xo(:,keep);
    cfi.yo=cfi.yo(:,keep);
    cf(i)=cfi;
end
end
