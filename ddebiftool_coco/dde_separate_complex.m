function [lam,trivlam,nontriv]=dde_separate_complex(lambda,triv)
nl=length(lambda);
ntriv=length(triv);
dist=abs(repmat(lambda(:).',ntriv,1)-repmat(triv(:),1,nl));
sel=munkres(dist);
sel=sel(sel>0);
nontriv=1:nl;
nontriv(sel)=[];
trivlam=lambda(sel);
lam=lambda;
lam(sel)=[];
end
