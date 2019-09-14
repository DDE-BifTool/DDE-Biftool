%% test high-order derivative
clear
freq=1;
f=@(x)real((exp(1i*freq*x)+exp(-1i*freq*x))/2);
order=3;
df=@(x)freq^order*0.5*real(1i^order*exp(1i*freq*x)+(-1i)^order*exp(-1i*freq*x));
gr=(sqrt(5)-1)/2;
ratio=2.^(1/order);
deriv=nmfm_Dbary(order,ratio,'n',order);
h=1./2.^(0:12);
xrg=linspace(-pi,pi,20);
dfa=NaN(length(h),length(xrg));
for i=1:length(h)
    for k=1:length(xrg)
        xdev=xrg(k)+h(i)*deriv.x(:)';
        fdev=f(xdev).';
        f0=f(0);
        dfa(i,k)=(deriv.D'*fdev+deriv.D0*f0)/h(i)^order-df(xrg(k));
    end
end
err_ratio=max(abs(dfa(1:end-1,:))./dfa(2:end,:),[],2)';
h_ratio=h(1:end-1)./h(2:end);
clf
semilogx(h(1:end-1),log10(err_ratio)./log10(h_ratio),'.-');
grid on
%% determine order of formula
deriv_orders=7;
polyorders=15;
D_orders=zeros(1,deriv_orders);
rng(4);
dfp=NaN(deriv_orders,polyorders);
jump=1e8;
hdev=1;
for i=1:deriv_orders
    ratio=3^(1/i);
    deriv=nmfm_Dbary(i,ratio,'n',i,'h',hdev);
    for k=1:polyorders
        p=rand(1,k+1);
        xe=rand(1);
        dp=p;
        for l=1:i
            dp=polyder(dp);
        end
        dfp(i,k)=(polyval(p,xe+deriv.x(:)'*hdev)*deriv.D+polyval(p,xe)*deriv.D0)/hdev^i-polyval(dp,xe);
        if k>1 && abs(dfp(i,k-1))~=0 && abs(dfp(i,k))>jump*abs(dfp(i,k-1))
            D_orders(i)=k-1;
            break
        end
    end
end
