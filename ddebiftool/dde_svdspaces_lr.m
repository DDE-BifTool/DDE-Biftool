function [v,w,evs]=dde_svdspaces_lr(A,nulldim,varargin)
%% determine left and right null space of large sparse matrix (possibly non-square)
% This routine relies on the property of the LU factorization that it puts
% the zero elements at the bottom of U
%
% $Id: dde_svdspaces_lr.m 309 2018-10-28 19:02:42Z jansieber $
%% 
default={'nulltol',1e-8};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[v,w]=dde_nullspaces_lr(A,'nulltol',options.nulltol,pass_on{:});
tolnulldim=size(v,2);
if tolnulldim==nulldim
    evs=[];
    return
end
evdim=nulldim-tolnulldim;
J1=[A,w;v',sparse(tolnulldim,size(w,2))];
[s1,s2]=size(J1);
Je=[sparse(s1,s1),J1; J1',sparse(s2,s2)];
[L,U,P,Q,R]=lu(Je);
%[Lt,Ut,Pt,Qt,Rt]=lu(J1');
opts.issym=1;
opts.isreal=1;
[V,D]=eigs(@(v)Q*(U\(L\(P*(R\v)))),...
    s1+s2,2*evdim,'sm',opts);
%[W,~]=eigs(@(v)(Qt*(Ut\(Lt\(Pt*(Rt\(Q*(U\(L\(P*(R\v)))))))))),...
%    size(J1,1),evdim,-1e-6,opts);
sv=diag(D);
[~,is]=sort(abs(sv));
sv=sv(is);
we=V(1:s1,is);
ve=V(s1+1:end,is);
sel=sv>=0;
evs=sv(sel);
ve=ve(1:size(A,2),sel);
we=we(1:size(A,1),sel);
[v,~]=qr([v,ve],0);%qr([v,V(1:size(A,2),:)],0);
[w,~]=qr([w,we],0);%qr([w,W(1:size(A,1),:)],0);
end
