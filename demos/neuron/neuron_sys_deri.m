%% User-provided partial derivatives of right-hand side |f|
%
% Implementation for tutorial demo neuron:
%
%   function J=neuron_sys_deri(xx,par,nx,np,v)
% 
% <html>
% $Id: neuron_sys_deri.m 367 2019-07-14 21:56:46Z jansieber $
% </html>
%
% If a user-provided function is not provided, DDE-Biftool will use
% finite-difference approximation (implemented in |df_deriv.m|).
%%
function J=neuron_sys_deri(xx,par,nx,np,v)

%% Parameter vector
% |par| $=[\kappa, \beta, a_{12}, a_{21},\tau_1,\tau_2, \tau_s]$.

J=[];
if length(nx)==1 && isempty(np) && isempty(v)
    %% First-order derivatives wrt to state nx+1
    if nx==0 % derivative wrt x(t)
        J(1,1)=-par(1);
        J(2,2)=-par(1);
    elseif nx==1 % derivative wrt x(t-tau1)
        J(2,1)=par(4)*(1-tanh(xx(1,2))^2);
        J(2,2)=0;
    elseif nx==2 % derivative wrt x(t-tau2)
        J(1,2)=par(3)*(1-tanh(xx(2,3))^2);
        J(2,2)=0;
    elseif nx==3 % derivative wrt x(t-tau_s)
        J(1,1)=par(2)*(1-tanh(xx(1,4))^2);
        J(2,2)=par(2)*(1-tanh(xx(2,4))^2);
    end
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% First order derivatives wrt parameters
    if np==1 % derivative wrt kappa
        J(1,1)=-xx(1,1);
        J(2,1)=-xx(2,1);
    elseif np==2 % derivative wrt beta
        J(1,1)=tanh(xx(1,4));
        J(2,1)=tanh(xx(2,4));
    elseif np==3 % derivative wrt a12
        J(1,1)=tanh(xx(2,3));
        J(2,1)=0;
    elseif np==4 % derivative wrt a21
        J(2,1)=tanh(xx(1,2));
    elseif np==5 || np==6 || np==7 % derivative wrt tau
        J=zeros(2,1);
    end;
elseif length(nx)==1 && length(np)==1 && isempty(v)
    %% Mixed state, parameter derivatives
    if nx==0 % derivative wrt x(t)
        if np==1 % derivative wrt beta
            J(1,1)=-1;
            J(2,2)=-1;
        else
            J=zeros(2);
        end;
    elseif nx==1 % derivative wrt x(t-tau1)
        if np==4 % derivative wrt a21
            J(2,1)=1-tanh(xx(1,2))^2;
            J(2,2)=0;
        else
            J=zeros(2);
        end;
    elseif nx==2 % derivative wrt x(t-tau2)
        if np==3 % derivative wrt a12
            J(1,2)=1-tanh(xx(2,3))^2;
            J(2,2)=0;
        else
            J=zeros(2);
        end;
    elseif nx==3 % derivative wrt x(t-tau_s)
        if np==2 % derivative wrt beta
            J(1,1)=1-tanh(xx(1,4))^2;
            J(2,2)=1-tanh(xx(2,4))^2;
        else
            J=zeros(2);
        end;
    end;
elseif length(nx)==2 && isempty(np) && ~isempty(v)
    %% Second order derivatives wrt state variables
    if nx(1)==0 % first derivative wrt x(t)
        J=zeros(2);
    elseif nx(1)==1 % first derivative wrt x(t-tau1)
        if nx(2)==1
            th=tanh(xx(1,2));
            J(2,1)=-2*par(4)*th*(1-th*th)*v(1);
            J(2,2)=0;
        else
            J=zeros(2);
        end;
    elseif nx(1)==2 % derivative wrt x(t-tau2)
        if nx(2)==2
            th=tanh(xx(2,3));
            J(1,2)=-2*par(3)*th*(1-th*th)*v(2);
            J(2,2)=0;
        else
            J=zeros(2);
        end
    elseif nx(1)==3 % derivative wrt x(t-tau_s)
        if nx(2)==3
            th1=tanh(xx(1,4));
            J(1,1)=-2*par(2)*th1*(1-th1*th1)*v(1);
            th2=tanh(xx(2,4));
            J(2,2)=-2*par(2)*th2*(1-th2*th2)*v(2);
        else
            J=zeros(2);
        end
    end
end
%% Otherwise raise error
% Raise error if the requested derivative does not exist
if isempty(J)
    error(['SYS_DERI: requested derivative nx=%d, np=%d, size(v)=%d',...
        'could not be computed!'],nx,np,size(v));
end
end
