%% test nfm_dfdx_scalar
clear
% first test scalar functions
freq=2;
f0=@(x)real((exp(1i*freq*x)+exp(-1i*freq*x))/2);
xshift=1;
f=@(x)[f0(x);f0(x+xshift)];
df0=@(x,order)freq^order*0.5*real(1i^order*exp(1i*freq*x)+(-1i)^order*exp(-1i*freq*x));
df=@(x,order)[df0(x,order);df0(x+xshift,order)];
xrg=linspace(-pi,pi,20);
mxorder=5;
for i=1:mxorder
    h(i,1)=1e-4;
    for k=1:length(xrg)
        h0=h(i,max(1,k-1));
        [dfa(:,i,k),dfaest(:,i,k),h(i,k)]=nmfm_dfdx_scalar(@(x)f(xrg(k)+x),i,'h',h0,...
            'isvectorized',true);
        err(i,k)=norm(dfa(:,i,k)-df(xrg(k),i),inf);
        errest(i,k)=norm(dfaest(:,i,k)-df(xrg(k),i),inf);
    end
end
err_est=squeeze(max(abs(dfa-dfaest),[],1));