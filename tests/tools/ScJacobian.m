%% Jacobian of f in x (2nd order finite difference)
%%
function df=ScJacobian(f,x,varargin)
default={'h',1e-5};
options=dde_set_options(default,varargin,'pass_on');
nc=length(x);
nr=length(f(x));
df=zeros(nr,nc);
for i=1:nc
    xu=x;
    xu(i)=xu(i)+options.h;
    xl=x;
    xl(i)=xl(i)-options.h;
    df(:,i)=(f(xu)-f(xl))/(2*options.h);
end
end
