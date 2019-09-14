function y=dirderi(i,x,p,vinp)
sel=@(x0,y0)bsxfun(@(x,y)x(y),x0,y0);
v01=@(t)sel(vinp(0,t),1);
v11=@(t)sel(vinp(1,t),1);
v21=@(t)sel(vinp(2,t),1);
v31=@(t)sel(vinp(3,t),1);
v41=@(t)sel(vinp(4,t),1);
v51=@(t)sel(vinp(5,t),1);
switch i
case 1
t1 = p(6) * x(1);
y = [-p(1) * v01(-t1 - p(3)) - p(2) * v01(-t1 - p(4)) - p(5) * v01(0)];
case 2
t1 = (p(6) * x(1));
y = [2 * p(6) * v01(0) * (p(1) * v11(-t1 - p(3)) + p(2) * v11(-t1 - p(4)))];
case 3
t1 = (p(6) * x(1));
t2 = v01(0);
y = [-3 * p(6) ^ 2 * t2 ^ 2 * (p(1) * v21(-t1 - p(3)) + p(2) * v21(-t1 - p(4)))];
case 4
t1 = (p(6) * x(1));
t2 = v01(0);
t3 = (t2 ^ 2);
y = [4 * p(6) ^ 3 * t2 * t3 * (p(1) * v31(-t1 - p(3)) + p(2) * v31(-t1 - p(4)))];
case 5
t1 = (p(6) * x(1));
t2 = v01(0);
t2 = t2 ^ 2;
t3 = (p(6) ^ 2);
y = [-5 * t3 ^ 2 * t2 ^ 2 * (p(1) * v41(-t1 - p(3)) + p(2) * v41(-t1 - p(4)))];
end
end
