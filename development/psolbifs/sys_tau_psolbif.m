function tau=sys_tau_psolbif(k,ntau1,sys_tau,x,p,varargin)
if k<ntau1
    tau=sys_tau(k,x,p,varargin{:});
elseif k==ntau1
    if isempty(varargin)
        tau=zeros([1,size(x,3)]);
    else
        nx=varargin{1};
        np=varargin{2};
        if length(nx)==1 && isempty(np)
            tau=zeros([1,xsize(1),nvec]);
        elseif length(nx)==2
            tau=zeros([xsize(1),xsize(1),nvec]);
        elseif  length(nx)==1 && length(np)==1
            tau=zeros([1,xsize(1),nvec]);
        end
    end
else
    tau=sys_tau(k-ntau1,varargin{:});
end
end
