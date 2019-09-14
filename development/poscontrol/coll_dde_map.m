function [coll,extmesh]=coll_dde_map(funcs,pt,varargin)
%% Map values on mesh in point structure to collocation points
%
%% Inputs:
% * funcs: problem definition  (sys_tau, sys_ntau, tpdel needed)
% * pt: point structure (psol or hcli) with collocation mesh, degree and
% profile, and period
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
% $Id$
%%
default={'wrapJ',true,'c_is_tvals',false,'c',[],'nderivs',1};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% Determine where residuals are computed: mesh tc 
% only compute res & jac at specified points?
if ~options.c_is_tvals
    % compute at all collocation points
    tcoarse=pt.mesh(1:pt.degree:end);    % boundaries of collocation intervals
    neqs=pt.degree*(length(tcoarse)-1);  % number of equations from collocation points
    if isempty(options.c)
        options.c=reshape(poly_gau(pt.degree),1,pt.degree);
    end
    h_int=diff(tcoarse); % lengths of collocation intervals
    tc=tcoarse(ones(length(options.c),1),1:end-1)+...
        options.c(ones(length(h_int),1),:)'.*...
        h_int(ones(length(options.c),1),:);
else
    neqs=length(options.c); % only evaluate res,Jac at times c
    tc=options.c;
end
tc=tc(:)';
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
    n=size(pt.profile,1);
    yarr=NaN(n,d,neqs);
end
%% generate W, W', W^(k-1) for each delay
% W{t_i}) is W for delay t_i-1 (t_i=1 corresponds to tau=0), W{t_i}{2}
% is W' at delay t_i-1.
W=cell(d,options.nderivs+1);
y=cell(d,options.nderivs+1);
pt_eva=str2func([pt.kind,'_eva']);
oneqs=ones(1,neqs);
par_neqs=pt.parameter(1,:,oneqs);
for k=0:options.nderivs
    for t_i=1:d
        [y{t_i,k+1},W{t_i,k+1}]=pt_eva(pt,tc-tau(:,t_i)'/pt.period,...
            'diff',k,'wrapJ',options.wrapJ,pass_on{:});
        if funcs.tp_del && t_i<d && k==0
            yarr(:,t_i,:)=reshape(y{t_i,1},n,1,neqs);
            tau(:,t_i+1)=funcs.sys_tau(t_i,yarr(:,1:t_i,:),par_neqs);
        end
    end
    if ~options.wrapJ
        [W(:,k+1),extmesh]=dde_psol_sparse_unwrap(W(:,k+1),pt.mesh,pt.degree);
    end
end
if options.wrapJ
    extmesh=pt.mesh;
end
coll=struct('mesh',tc,'profile',{y},'jac',{W},'tau',tau);
end
