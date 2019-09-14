function J=rot_deriv(funcs,xx,par,nx,np,v)
%% adapt user provided derivative to rotating coordinates
%% UNTESTED!! Does not yet work.
%
% $Id: rot_deriv.m 309 2018-10-28 19:02:42Z jansieber $
%
A=funcs.rotation;
expA=funcs.exp_rotation;
omega=par(end);
ind_omega=length(par);
userpar=par(1:end-1);
xxrot=xx;
tau_ind=funcs.sys_tau();
tau=[0,userpar(tau_ind)];
for i=2:size(xx,2)
    xxrot(:,i)=expA(-omega*p(tau_ind(i)))*xxrot(:,i);
end
orig_funcs=funcs;
orig_funcs.sys_rhs=funcs.orig_rhs;
J=[];
if length(nx)==1 && isempty(np) && isempty(v),
    % first order derivatives of the state:
    J=funcs.orig_deri(orig_funcs,xxrot,userpar,nx,np,v);
    J=J*expA(-omega*tau(nx+1));
    if nx==0
        J=J-A*omega;
    end
elseif isempty(nx) && length(np)==1 && isempty(v),
    % first order parameter derivatives:
    if np==ind_omega
        % derivative wrt omega
        J=-A*xxrot(:,1);
        for i=2:size(xx,2)
            Ji=funcs.orig_deri(orig_funcs,xxrot,userpar,i-1,[],[]);
            Ji=Ji*expA(-omega*tau(i))*(-A*tau(i))*xxrot(:,i);
            J=J+Ji;
        end
    else
        J=funcs.orig_deri(orig_funcs,xxrot,userpar,nx,np,v);
        % change derivatives if p(np) is a delay
        tausel=find(ismember(np,tau_ind),1);
        if ~isempty(tausel)
            Jxtau=funcs.orig_deri(orig_funcs,xxrot,userpar,tausel,[],[]);
            J=J+Jxtau*expA(-omega*tau(tausel+1))*(-A*omega)*xxrot(:,i);
        end
    end
elseif length(nx)==2 && isempty(np) && ~isempty(v),
    % second order state derivatives
    nx=sort(nx);
    J=funcs.orig_deri(orig_funcs,xxrot,userpar,nx,np,v);
    J=J*expA(-omega*(tau(nx(1)+1)+tau(nx(2)+1)));
elseif length(nx)==1 && length(np)==1 && isempty(v),
    % mixed state parameter derivatives
    if np==ind_omega
        if nx==0 % deriv wrt to instant state
            J=-A;
            for i=2:size(xxrot,2)
                v=expA(-omega*tau(i))*(-A*tau(i))*xxrot(:,i);
                Ji=funcs.orig_deri(orig_funcs,xxrot,userpar,[0,i-1],[],v);
                J=J+Ji;
            end
        else % deriv wrt to delayed state
            J=funcs.orig_deri(orig_funcs,xxrot,userpar,nx,[],[]);
            J=J*expA(-omega*tau(nx+1))*(-A*tau(nx+1));
            for i=2:size(xxrot,2)
                v=expA(-omega*tau(i))*(-A*tau(i))*xxrot(:,i);
                Ji=funcs.orig_deri(orig_funcs,xxrot,userpar,[nx,i-1],[],v);
                Ji=Ji*expA(-omega*tau(nx+1));
                J=J+Ji;
            end
        end
    else
        J=funcs.orig_deri(orig_funcs,xxrot,userpar,nx,np,v);
        J=J*expA(-omega*tau(nx+1));
        % change derivatives if p(np) is a delay
        tausel=find(ismember(np,tau_ind),1);        
        if ~isempty(tausel)
            v=expA(-omega*tau(tausel+1))*(-A*omega)*xxrot(:,tausel+1);
            Jxxtau=funcs.orig_deri(orig_funcs,xxrot,userpar,[nx,tausel],[],v);
            Jxxtau=Jxxtau*expA(-omega*tau(nx+1));
            J=J+Jxxtau;
            if tausel==nx+1
                Jtau=funcs.orig_deri(orig_funcs,xxrot,userpar,nx,[],[]);
                Jtau=Jtau*expA(-omega*tau(tausel+1))*(-A*omega);
                J=J+Jtau;
            end
        end
    end
end

if isempty(J)
    error('SYS_DERI: requested derivative (nx,np,v)=(%d,%d,%d) does not exist!',...
        nx,np,size(v,1));
end

end
