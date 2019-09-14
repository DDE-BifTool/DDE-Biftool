function cf=dde_chainrule_combinatorics(maxorder)
%% differentiate g(x(s),y(s)) 1..maxorder (1..p) times
% Results is sum of i*D[j,k](g){D[1]x^(m(1),...D[p]x^m(p), D[1]y^(n(1),...D[p]y^n(p)}
% coefficient for qth derivativeis stored at
% i=cf(q).fac
% j=cf(q).gx;
% m(r)=cf(q).xo;
cf=repmat(struct('fac',1,'gx',1,'xo',[1;zeros(maxorder-1,1)]),maxorder,1);
for i=1:maxorder-1
    nc=cf(i);
    numc=length(nc.fac);
    nrnew=sum(nc.xo(:)>0)+numc;
    xo=NaN(maxorder,nrnew);
    fac=NaN(1,nrnew);
    gx=NaN(1,nrnew);
    ind=0;
    for j=1:numc
        %% differentiate g wrt x, add 1st deriv of x
        ind=ind+1;
        fac(ind)=nc.fac(j);
        gx(ind)=nc.gx(j)+1;
        xo(:,ind)=nc.xo(:,j);
        xo(1,ind)=xo(1,ind)+1;
        ix=find(nc.xo(:,j)>0);
        for k=1:length(ix)
            ind=ind+1;
            fac(ind)=nc.fac(j);
            gx(ind)=nc.gx(j);
            xo(:,ind)=nc.xo(:,j);
            %% differentiate any derivative of x once more
            xo(ix(k)+1,ind)=xo(ix(k)+1,ind)+1;
            fac(ind)=fac(ind)*xo(ix(k),ind);
            xo(ix(k),ind)=xo(ix(k),ind)-1;
        end
    end
    sxo=2;
    %% check for equal entries
    mat=[gx;xo;fac]';
    mat=sortrows(mat);
    [~,iu]=unique(mat(:,1:end-1),'rows');
    iu(end+1)=size(mat,1)+1; %#ok<AGROW>
    nfac=zeros(1,length(iu)-1);
    for j=1:length(iu)-1
        nfac(j)=sum(mat(iu(j):iu(j+1)-1,end));
    end
    mat=[mat(iu(1:end-1),1:end-1),nfac'];
    mat=sortrows(mat,[1,maxorder-1:1]);
    mat=mat';
    cf(i+1).fac=mat(end,:);
    cf(i+1).gx=mat(1,:);
    cf(i+1).xo=mat(sxo+(0:maxorder-1),:);
end
for i=1:maxorder
    if isfield(cf(i),'xo')
        cf(i).xo=cf(i).xo(1:i,:);
    end
end
end
