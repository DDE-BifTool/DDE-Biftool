function [x,converged,J]=NSolve(f,x0,varargin)
defaults={'maxit',10,'tolerance',1e-5,...
    'print',1,'damping',1,'damping_threshold',1.2,...
    'update_threshold',2e-1,'discard_threshold',1e1,'maxcor',Inf,...
    'cnorm',@(x)norm(x,inf),'numjac',false};
[options,pass_on]=dde_set_options(defaults,varargin,'pass_on');
x=x0;
converged=0;
oldcornorm=options.maxcor;
Jup2date=false;
for i=1:options.maxit
    if ~Jup2date
        if ~options.numjac
            [y,J]=f(x);
        else
            y=f(x);
            J=MyJacobian(f,x,'f0',y,pass_on{:});
        end
        Jfresh=true;
    else
        y=f(x);
        Jfresh=false;
    end
    cor=-J\y;
    ynorm=options.cnorm(y);
    cornorm=options.cnorm(cor);
    if isnan(cornorm)
            converged=0;
            break
    end
    if oldcornorm*options.update_threshold>cornorm && i>1
        % update Jacobian if convergence is poor
        Jup2date=true;
    else
        Jup2date=false;
    end
    if i<=1 || cornorm<oldcornorm*options.discard_threshold
        % discard correction if convergence is really poor
        x=x+options.damping*cor;
        if options.damping<1 &&  i>1 && ...
                cornorm<oldcornorm*options.damping*options.damping_threshold
            options.damping=(1+options.damping)/2;
        end
        oldcornorm=cornorm;
        acc=true;
    else
        acc=false;
        if Jfresh
            converged=0;
            break
        end
    end
    if options.print>0
        fprintf('it=%d, |cor|=%g, |res|=%g, acc=%d, renew J=%d\n',...
            i,cornorm,ynorm,acc,~Jup2date);
    end
    if max(cornorm,ynorm)<options.tolerance
        converged=1;
        break
    end
    if max(cornorm,ynorm)>options.maxcor
        converged=0;
        break
    end
end
% if nargout>2 && converged
%     [y,J]=f(x); %#ok
% end
if converged
    converged=i;
end
end
%% Jacobian of f in x (2nd order finite difference)
%%
function df=MyJacobian(f,x,varargin)
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
