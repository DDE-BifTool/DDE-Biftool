%% Jacobian of f in x (2nd order finite difference)
%%
function df=ScJacobian(f,x,h)
nc=length(x);
nr=length(f(x));
df=zeros(nr,nc);
for i=1:nc
    xu=x;
    xu(i)=xu(i)+h;
    xl=x;
    xl(i)=xl(i)-h;
    df(:,i)=(f(xu)-f(xl))/(2*h);
end
end
