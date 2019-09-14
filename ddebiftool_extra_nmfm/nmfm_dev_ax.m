function fun=nmfm_dev_ax(avec,funs)
%% linear combination of history functions
%
% $Id: nmfm_dev_ax.m 309 2018-10-28 19:02:42Z jansieber $
%%
nf=length(avec);
for i=1:nf
    funs(i).v=avec(i)*funs(i).v;
end
vvec=cat(2,funs.v);
lvec=cat(2,funs.lambda);
tvec=cat(2,funs.t);
fun=nmfm_dev_fun(vvec,'lambda',lvec,'t',tvec);
end