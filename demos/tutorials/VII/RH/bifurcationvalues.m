function [xstar,Iapp,tau] = bifurcationvalues(a,b,c,d,chi,r,S)
    xstar=1/(3*a)*(b-d+sqrt((b-d)^2-3*a*c*S));
    Iapp=xstar^2*(a*xstar-b+d)+c*(S*(xstar-chi)-1);
    % auxiliary variables to calculate tau
    A=xstar^2*((3*a*xstar+2*d)^2-4*b^2)-2*r*xstar*(2*d*xstar-1)*(3*a*xstar-2*b+2*d)+...
        r^2*(4*d*xstar*(-2*b*xstar+d*xstar-1)+1);
    B=9*a^2*xstar^4+2*r*xstar*(3*a*xstar-2*b+2*d)-4*b^2*xstar^2-4*d*xstar+r^2+1;
    omega=sqrt(-1/2*B-1/2*sqrt(B^2-4*A));
    if ~isreal(omega)
        omega=sqrt(-1/2*B+1/2*sqrt(B^2-4*A));
    end
    Y=omega/(2*b)*(-1/xstar+2*d/(1+omega^2)+r*(2*b-2*d-3*a*xstar)/(r^2+omega^2));
    tau=1/omega*(asin(Y)+2*pi);
end