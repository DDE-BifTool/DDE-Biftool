function [r,J]=sys_cond_POBP(p,M)
[ncond,dim]=size(M);
J=repmat(p_axpy(0,p,[]),ncond,1);
ot=ones(length(p.mesh),1);
for i=length(J):-1:1
    p2=p_axpy(0,p,[]);
    p2.profile=0*p2.profile;
    p2.profile(1:dim,:)=M(i*ot,:).';
    [r(i),J(i)]=p_dot(p,p2);
end
end
