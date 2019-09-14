function err=struct_error(s1,s2)
tol=eps;
tolmax=1e0;
d1=1;
while ~isempty(d1) && tol<tolmax
    [~,d1]=comp_struct(s1,s2,1,0,tol);
    tol=2*tol;
end
err=tol;
end