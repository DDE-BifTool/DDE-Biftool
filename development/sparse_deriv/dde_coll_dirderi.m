function varargout=dde_coll_dirderi(funcs,pt,free_par,varargin)
default={'Dtmat',[],'order',1,'nf',size(pt.profile,1),...
    'coll',[],'xdev',[],'output','res+J',...
    'pattern',{@dde_jacobian1pattern,@dde_jacobian2pattern}};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
ind=dde_ind_from_point(pt,free_par);
%% Dtmat is matrix in front of leading derivative
if isempty(options.Dtmat)
    options.Dtmat=eye(options.nf);
end
%% get values at collocation points
if isempty(options.coll)
    [coll,extmesh]=dde_coll_map(funcs,pt,'nderivs',options.order,...
        'free_par_ind',free_par,pass_on{:});
end
%% matrix dimensions
dim=size(pt.profile,1);
nc=size(coll.jac{1},1)/dim; % num. coll. points
d=size(coll.jac,1);   % number of delays +1
%% base point
u=dde_x_from_point(pt,free_par);
nu=length(u);
x=u(ind.profile(:));
p=pt.parameter;
T=u(ind.period);
%% r.h.s. and derivatives
E0=cat(1,coll.jac{:,1});
E10=coll.jac{1,2}/T;
xc=E0*x;
dfxp_fun=@(ord,dy,dq)dde_coll_sysvaldir(ord,funcs,xc,p,...
    dy,q_ins(dq,0*p,free_par),[dim,nc,d]);
%% determine residual if required
res=[];
Dts=sparse_blkdiag(options.Dtmat(:,:,ones(1,nc)));
if any(strcmp(options.output,{'res+J','res','J+res'}))
    F=dfxp_fun(0,[],[]);
    res=Dts*E10*x-F;
    if strcmp(options.output,'res')
        varargout=output(options.output,res,[],[],[]);
        return
    end
end
%% determine pattern of 1st order Jacobian
[ir,ic]=find(E0);
E0shape=sparse(ir,ic,ones(size(ir)),size(E0,1),size(E0,2));
dfshape=repmat(kron(speye(nc),ones(options.nf)),1,d)*E0shape;
dfshape=cat(2,dfshape,ones(size(dfshape,1),1+length(free_par)));
%% deviations needed to determine Jacobian
Jpat=options.pattern{options.order}(dfshape,'ncol',nu);
xdev=Jpat.dev;
y=xdev(ind.profile(:),:);
q=xdev(ind.parameter,:);
S=xdev(ind.period,:)/T;
%% dependence of tau on free parameters
itau=funcs.sys_tau();
tau=[0;pt.parameter(itau)'];
[istau,indpt]=ismember(free_par,itau);
dtaudp=full(sparse(indpt(istau)+1,find(istau),ones(sum(istau),1),...
    itau+1,length(free_par)));
%% expand delays to constant functions
to_c=@(th)reshape(repmat(reshape(th,1,[]),dim*nc,1),dim*nc*d,[]);
%% First order deviation in x
E1=cat(1,coll.jac{:,2})/T;
bdiag=@(vec)sparse_blkdiag(reshape(vec,dim*nc,1,[]));
dtdev=dtaudp*q-tau*S;
dY1=E0*y-bdiag(E1*x)*dtdev;
dfxp=dfxp_fun(options.order,dY1,q);
if options.order==1
    Jvec=Dts*E10*y-Dts*E10*x*S-dfxp;
    J=Jpat.expand(Jvec);
    varargout=output(options.output,res,J,extmesh,coll);
    return
end
%% second-order Jacobian requested
S_d=S(ones(d,1),:);
S_n=S(ones(nc*dim,1),:);
E2=cat(1,coll.jac{:,3})/T^2;
dt2dev=tau*(S.^2)-(dtaudp*q).*S_d;
dY2=-2*(E1*y).*to_c(dtdev)+bdiag(E2*x)*(dtdev.^2)-2*bdiag(E1*x)*dt2dev;
dfxp1=dfxp_fun(1,dY2,0*q);
Jvec=-2*(Dts*E10*y).*S_n+2*(Dts*E10*x).*(S_n.^2)-dfxp-dfxp1;
J=Jpat.expand(Jvec);
varargout=output(options.output,res,J,extmesh,coll);
end
%%
function qext=q_ins(q,qframe,free_par)
ndev=size(q,2);
qext=zeros(1,length(qframe),ndev);
qext(1,free_par,:)=reshape(full(q),1,length(free_par),ndev);
end
%% output
function out=output(opt,res,J,extmesh,coll)
switch opt
    case 'res'
        out={res};
    case {'J','jac','Jacobian'}
        out={J,extmesh,coll};
    case {'res+J','res+jac','res+Jacobian'}
        out={res,J,extmesh,coll};
    case {'J+res'}
        out={J,res,extmesh,coll};
end
end