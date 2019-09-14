function ind=dde_coll_subinterval(msh,t,orient)
%% return ind s.t. msh(ind(i))<=t(i)<msh(ind(i)+1) (if orient==0, default)
% msh(ind(i))<t(i)<=msh(ind(i)+1) (if orient==1)
%
% convention ind(i)=0 if t<msh(1), ind(i)=length(msh)+1 if t>msh(end).
% If t(i)=msh(end) then ind(i)=length(msh) for orient=0, and, if
% t(i)=msh(1) then ind(i)=1 for orient=1.
%
% msh is assumed to be sorted, t may be not (but will be sorted and
% reverse sorted inside if not)
%
% consecutive entries in msh are assumed to be different
%
% $Id: dde_coll_subinterval.m 369 2019-08-27 00:07:02Z jansieber $
%%
if nargin<3
    orient=0;
end
nt=length(t);
nm=length(msh);
irev=1:nt;
if ~issorted(t)
    [t,irev]=sort(t);
end
ind=zeros(size(t));
jt=find(t<msh(1),1,'last');
if isempty(jt)
    jt=1;
end
while jt<=nt && t(jt)==msh(1)
    ind(jt)=1;
    jt=jt+1;
end
jm=2;
while jt<=nt && jm<=nm
    if (~orient && t(jt)<msh(jm)) || (orient &&  t(jt)<=msh(jm))
        ind(jt)=jm-1;
        jt=jt+1;
    else
        jm=jm+1;
    end
end
while jt<=nt && t(jt)==msh(nm)
    ind(jt)=nm;
    jt=jt+1;
end
ind(jt:end)=nm+1;
ind(irev)=ind;
end
