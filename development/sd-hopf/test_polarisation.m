%% test polarization identities
clear
order=2;
ntests=10;
rng(3);
zr0=randn(1,ntests);
zi0=randn(1,ntests);
[m,f]=nmfm_pol_real_from_complex(order);
for i=1:ntests
    combs0(i,:)=(zr0(i)*m(1,:)+zi0(i)*m(2,:)).^order.*f;
    scomb0(i)=sum(combs0(i,:));
    pr0(i)=(zr0(i)+1i*zi0(i))^order;
end
assert(all(abs(scomb0-pr0)<1e-10));
%% general polarization, no grouping
zr=rand(order,ntests);
zi=randn(order,ntests);
isreal=randi(3,order,ntests)<=1;
zi(isreal)=0;
for i=1:ntests
    dev=(zr(:,i)+1i*zi(:,i)).';
    [x{i},fac{i}]=nmfm_pol_from_devs(dev,order,num2cell(1:order));
    comb{i}=x{i}.^order.*fac{i};
    scomb(i)=sum(comb{i});
    pr(i)=prod(dev);
end
[maxerr,imxerr]=max(abs(pr-scomb));
maxlen=max(cellfun(@(x)length(x),comb));
fprintf('max length=%d, pr(%d)-scomb(%d)=%g\n',maxlen,imxerr,imxerr,maxerr);
%% general polarization, random grouping
zr=rand(order,ntests);
zi=randn(order,ntests);
isreal=randi(3,order,ntests)<=1;
zi(isreal)=0;
ngroups=randi(order,1,ntests);
irg=1:order;
for i=ntests:-1:1
    seteq=[0,sort(randperm(order-1,ngroups(i)-1)),order]+0.5;
    for k=1:length(seteq)-1
        ind=irg(irg>seteq(k)&irg<seteq(k+1));
        if isempty(ind)
            continue
        end
        zr(ind,i)=zr(ind(1),i);
        zi(ind,i)=zi(ind(1),i);
    end
    shuffle=randperm(order);
    zr(:,i)=zr(shuffle,i);
    zi(:,i)=zi(shuffle,i);
end
%%
for i=1:ntests
    dev=(zr(:,i)+1i*zi(:,i)).';
    groups{i}=nmfm_dev_group(dev);
    [x{i},fac{i}]=nmfm_pol_from_devs(dev,order,groups{i});
    comb{i}=x{i}.^order.*fac{i};
    scomb(i)=sum(comb{i});
    pr(i)=prod(dev);
end
[maxerr,imxerr]=max(abs(pr-scomb));
maxlen=max(cellfun(@(x)length(x),comb));
fprintf('max length=%d, pr(%d)-scomb(%d)=%g\n',maxlen,imxerr,imxerr,maxerr);