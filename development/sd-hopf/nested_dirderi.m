function y=dirderi(i,x,vinp,p)
sel=@(x0,y0)bsxfun(@(x,y)x(y),x0,y0);
v01=@(t)sel(vinp(0,t),1);
v11=@(t)sel(vinp(1,t),1);
v21=@(t)sel(vinp(2,t),1);
v31=@(t)sel(vinp(3,t),1);
v41=@(t)sel(vinp(4,t),1);
v51=@(t)sel(vinp(5,t),1);
switch i
case 1
y = [-v01(-x(1))];
case 2
y = [2 * v11(-x(1)) * v01(0)];
case 3
t1 = v01(0);
y = [-3 * v21(-x(1)) * t1 ^ 2];
case 4
t1 = v01(0);
t2 = (t1 ^ 2);
y = [4 * v31(-x(1)) * t1 * t2];
case 5
t1 = v01(0);
t1 = t1 ^ 2;
y = [-5 * v41(-x(1)) * t1 ^ 2];
end
end
