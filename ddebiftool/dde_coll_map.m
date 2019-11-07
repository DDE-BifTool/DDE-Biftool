function [coll,extmesh]=dde_coll_map(funcs,pt,varargin)
%% Map values on mesh in point structure to collocation points
%
%% Inputs:
% * funcs: problem definition  (sys_tau, sys_ntau, tpdel needed)
% * pt: point structure (eg, psol or hcli) with collocation mesh, degree and
% profile, and period
% * mesh evaluation function
%% Optional (name-value pairs):
% * wrapJ (true): wrap evaluation of points around periodically and adjust
% Jacobian correspondingly (for psol)
% * c_is_tvals: (false) whether collocation points are given explicitly in
% relation to the entire period
% * c (1 x pt.degree or 1 x neqs array): if c_is_tvals is false then tese
% are the points where collocation interpolation is evaluated. Otherwise,
% these are the collocation points in one subinterval (unscaled). If c is
% empty then the Gauss-Legendre points are used.
% * nderivs (default 1) number of derivatives to be computed
%
%% Outputs:
% * y (nderivs+1) x (ntau+1) cell array of n x neqs interpolation values
% * W (nderivs+1) x (ntau+1) cell array of (n*neqs) x numel(pt.profile)
% sparse Jacobian matrices
% * tc (1 x neqs) mesh of collocation points used (equal c if c is
% nonempty and c_is_tvals is true)
% * tau neqs x (ntau+1) array of delays
% * extmesh (if wrapJ is false, this contains the mesh points for all
%  all columns of y{k,d}.
%
% $Id: dde_coll_map.m 369 2019-08-27 00:07:02Z jansieber $
%%
default={'wrapJ',true,...
    'c_is_tvals',false,'c',[],'submesh_limit',0,...
    'nderivs',1,'deriv_order',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% Determine where residuals are computed: mesh tc 
% only compute res & jac at specified points?
n=size(pt.profile,1);
if ~options.c_is_tvals
    if isempty(options.c) || ischar(options.c)
        options.c=dde_coll_set_grid('collocation',pt.degree,'type',options.c,...
            'lhs_num',funcs.lhs_matrix(n));
    end
    if 1==max(options.c)
        options.submesh_limit=1;
    end
    tc=dde_coll_meshfill(pt.mesh(1:pt.degree:end),pt.degree,...
        'purpose','collocation','grid',options.c);
else
    tc=options.c(:)';
end
neqs=length(tc);
%% number of delays
% array tau will contain all delays at all collocation points
if ~funcs.tp_del
    n_tau=funcs.sys_tau();           % delay numbers
    tau=pt.parameter(n_tau);         % delay values
    d=length(n_tau)+1;               % number of delays
    tau=repmat([0;tau(:)]',neqs,1);  % # delays (incl tau=0) x evaluations
else
    n_tau=funcs.sys_ntau();
    d=n_tau+1;
    tau=[zeros(neqs,1),NaN(neqs,n_tau)];
    yarr=NaN(n,d,neqs);
end
if isempty(options.deriv_order)
    deriv_order=zeros(1,d);
else
    deriv_order=options.deriv_order;
end
%% generate W, W', W^(k-1) for each delay
% W{t_i}) is W for delay t_i-1 (t_i=1 corresponds to tau=0), W{t_i}{2}
% is W' at delay t_i-1.
W=cell(d,options.nderivs+1);
y=cell(d,options.nderivs+1);
oneqs=ones(1,neqs);
par_neqs=pt.parameter(1,:,oneqs);
for t_i=1:d
    for k=0:options.nderivs
        [y{t_i,k+1},W{t_i,k+1}]=pt_eval(pt,tc-tau(:,t_i)'/pt.period,'diff',k,...
            'wrapJ',options.wrapJ,'submesh_limit',options.submesh_limit,pass_on{:});
        if funcs.tp_del && t_i<d && k==deriv_order(t_i)
            yarr(:,t_i,:)=reshape(y{t_i,deriv_order(t_i)+1},n,1,neqs);
            tau(:,t_i+1)=funcs.sys_tau(t_i,yarr(:,1:t_i,:),par_neqs);
        end
    end
end
if ~options.wrapJ
    [W(:),extmesh]=pt_sparse_unwrap(pt,W(:));
end
if options.wrapJ
    extmesh=pt.mesh;
end
coll=struct('mesh',tc,'profile',{y},'jac',{W},'tau',tau);
end
%%
function [y,W]=pt_eval(pt,x,varargin)
% INPUT:
% * pt with profile: profile on mesh t, of degree order polynomials
% * x: point(s) where to evaluate
% OUTPUT:
%	px value of profile at x
%   J: sparse Jacobian:  y(:)=J*profile(:) if wrapJ is true, otherwise
%   J is structure with row indices, col indices and values for sparse
%   matrix, but which may contain negative column indices
default={'wrapJ',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% coarse mesh
%% wrap requested x into pt.mesh
bd=pt.mesh([1,end]);
len=bd(2)-bd(1);
pt.mesh=(pt.mesh-bd(1))/len;
x=(x-bd(1))/len;
c=x-floor(x);
c(x==1)=1; % for rotations & non-periodic collocation profiles
[y,W]=dde_coll_eva(pt.profile,pt.mesh,c,pt.degree,'kron',true,pass_on{:});
if ~options.wrapJ
    %% shift column indices toward past if x was negative
    n=size(pt.profile,1);
    nt=length(pt.mesh)-1;
    ishift=floor(x(:));
    sel1=x==floor(x)&x>0;
    ishift(sel1)=ishift(sel1)-1;
    ishiftn=reshape(repmat(ishift,1,n)',[],1);
    [ir,ic,Jvals]=find(W);
    %% determine column indices of points outside base interval
    ic=ic+ishiftn(ir)*n*nt;
    imin=min(ishift)*n*nt+1; % minimal index (whole base intervals)
    subint=floor((min(ic)-imin)/pt.degree);
    imin=imin+subint*pt.degree;
    %% store structure information for Jacobian
    W=struct('ir',ir,'ic',ic,'vals',Jvals,...
        'nr',size(W,1),'nc',size(W,2),'x',x*len+bd(1),'nkron',n);
end
end
%% convert structs to sparse matrices for psol
function [W,mesh]=pt_sparse_unwrap(pt,W)
% W{t_i} are structs, check for minimal index and shift to create
% sparse matrices
imin=min(cellfun(@(w)min(w.ic),W));
imax=max(cellfun(@(w)max(w.ic),W));
xmin=min(cellfun(@(w)min(w.x),W));
xmax=max(cellfun(@(w)max(w.x),W));
%% shift matrix indices and create sparse matrices
for t_i=1:length(W)
    W{t_i}.ic=W{t_i}.ic-imin+1;
    W{t_i}=sparse(W{t_i}.ir,W{t_i}.ic,W{t_i}.vals,...
        W{t_i}.nr,imax-imin+1);
end
%% extend mesh periodically, assuming pt.mesh([1,end])==[0,1]
assert(pt.mesh(1)==0&&pt.mesh(end)==1);
xfmin=floor(xmin);
xfmax=floor(xmax);
cmesh=pt.mesh(1:end);
% extend mesh by whole units of period to cover all x in W
cmesh=repmat(cmesh(:),1,xfmax-xfmin+1)+repmat(xfmin:xfmax,length(cmesh),1);
c_is_int=cmesh==round(cmesh);
cmesh(1,c_is_int(1,:)&cmesh(1,:)>=1)=NaN;
cmesh(end,c_is_int(end,:)&cmesh(end,:)<=0)=NaN;
cmesh=cmesh(:)';
cmesh=cmesh(~isnan(cmesh));
% restrict to necessary subintervals
icmin=find(cmesh<=xmin,1,'last');
icmin=floor((icmin-1)/pt.degree)*pt.degree+1;
icmax=find(cmesh>=xmax,1,'first');
icmax=ceil((icmax-1)/pt.degree)*pt.degree+1;
mesh=cmesh(icmin:icmax);
end
%%
