function df=dde_sparsejacobian(f,x,varargin)
%% compute Jacobian using finite differences (optionally prescribed pattern)
% optional inputs: 
%
% * 'h': deviation (default 1e-5)
% * 'order': order of finite difference formulas (1 or 2, default 2)
% * 'f0': f(x) (if f(x) has already been computed, default [] to recompute)
% * 'vectorised', if f can be called with several vectors (default false)
% * 'pattern': struct with fields 'color', 'column_bd', 'rows'
%
%% read in optional arguments
defaults={'h',1e-6,'order',2,'f0',[],'vectorized',false,'pattern',[],...
    'isdfx',false,'xrep',false};
options=dde_set_options(defaults,varargin,'pass_on');
if ~options.isdfx
    dfx=@(x,v)dfxloc(f,x,v,options.h,options.order,options.f0,options.vectorized);
else
    dfx=f;
end
nx=length(x);
pat=options.pattern;
%% which deviations to take
if isempty(pat)
    nf=size(dfx(x,0*x),1);
    nx=size(x,1);
    pat=dde_jacobianpattern(sparse(ones(nf,nx)));
end
ndev=max(pat.deviation);
%dev=full(sparse(pat.colors,pat.entries.columns(pat.column_bd(:,1)),ones(length(pat.colors),1)))';
dev=full(min(1,sparse(pat.i_sorted(:,2),pat.deviation,ones(length(pat.deviation),1))));
%% evaluate derivatives
if ~options.isdfx
    fdiffs=dfx(x,dev);
else
    if options.vectorized
        if options.xrep
            x=x(:,ones(1,ndev));
        end
        fdiffs=dfx(x,dev);
    else
        for i=ndev:-1:1
            fdiffs(:,i)=dfx(x,dev(:,i));
        end            
    end
end
nf=size(fdiffs,1);
%% insert values into sparse Jacobian
df=sparse(nf,nx);
ldf_ind=sub2ind(size(df),pat.i_sorted(:,1),pat.i_sorted(:,2));
ldev_ind=sub2ind([nf,ndev],pat.i_sorted(:,1),pat.deviation);
df(ldf_ind)=fdiffs(ldev_ind);
%% insert values into sparse Jacobian
% clr_bd=[1;find(diff(pat.colors))+1];
% clr_bd=[clr_bd,[clr_bd(2:end)-1;length(pat.colors)]];
% color_entries_bd=[pat.column_bd(clr_bd(:,1)),...
%     pat.column_bd(clr_bd(:,2),2)];
% dfvals=NaN(length(pat.entries.rows),1);
% dfind=0;
% for i=1:size(color_entries_bd,1)
%     ir=pat.entries.rows(color_entries_bd(i,1):color_entries_bd(i,2));
%     dfvals(dfind+(1:length(ir)))=fdiffs(ir,i);
%     dfind=dfind+length(ir);
% end
% df=sparse(pat.entries.rows,pat.entries.columns,dfvals);
end
%%
function fdiffs=dfxloc(f,x,dev,h,order,f0,vectorized)
ndev=size(dev,2);
dev=dev*h;
xu=x(:,ones(1,ndev))+dev;
if order==2
    xl=x(:,ones(1,ndev))-dev;
end
if vectorized
    if order==2
        fvals=f([xu,xl]);
        fdiffs=(fvals(:,1:ndev)-fvals(:,ndev+1:end))/(2*h);
    else
        fvals=f([xu,x]);
        fdiffs=(fvals(:,1:ndev)-fvals(:,ndev*ones(1,ndev)))/options.h;
    end
else
    if order==2
        for i=ndev:-1:1
            fdiffs(:,i)=(f(xu(:,i))-f(xl(:,i)))/(2*h);
        end
    else
        if isempty(f0)
            f0=f(x);
        end
        for i=ndev:-1:1
            fdiffs(:,i)=(f(xu(:,i))-f0)/h;
        end
    end
end
end
