function [y,W]=dde_coll_profile_dot(p1,p2,varargin)
%% integral int_a^b p1(s)'*p2(s) ds
% second output is matrix s.t. p1.profile(:)'*W*p2.profile(:)=y, optional
% arguments for other boundaries than mesh([1,end]), and selected
% components only
%
% $Id: dde_coll_profile_dot.m 369 2019-08-27 00:07:02Z jansieber $
%%
default={'derivatives',[0,0],'ind_comp',1:size(p1.profile,1),'bd',[0,1],...
    't_quad',[],'w_quad',[]};
options=dde_set_options(default,varargin,'pass_on');
%% treat case of int_a^a
y=0;
if options.bd(1)==options.bd(2)
    W=sparse(numel(p1.profile),numel(p2.profile));
    return
end
if isempty(options.t_quad)
    [options.t_quad,options.w_quad]=legpts(max(p1.degree,p2.degree),[0,1]);
end
%% store sign of integration direction
sgn=1-2*(options.bd(2)<options.bd(1));
bd=sort(options.bd);
%% join meshes & select subset between bd(1) and bd(2)
t1coarse=p1.mesh(1:p1.degree:end);
t2coarse=p2.mesh(1:p2.degree:end);
[tcoarse,~,ibd]=unique([bd,t1coarse,t2coarse]);
tcoarse=tcoarse(ibd(1):ibd(2));
%% create weights and nodes for integration (Gauss-Legendre)
% on each subinterval of joint mesh wmat is diagonal matrix of weights, t
% is sequence of nodes
w=dde_coll_meshfill(tcoarse,1,'grid',options.w_quad,'acc',false);
t=dde_coll_meshfill(tcoarse,1,'grid',options.t_quad);
nw=length(w);
dim=size(p1.profile,1);
dg=sparse(options.ind_comp,options.ind_comp,ones(length(options.ind_comp),1),dim,dim);
wmat=kron(sparse(1:nw,1:nw,w),dg);
[y1,W1]=dde_coll_eva(p1.profile,p1.mesh,t,p1.degree,...
    'diff',options.derivatives(1),'kron',true); %#ok<ASGLU>
[y2,W2]=dde_coll_eva(p2.profile,p2.mesh,t,p2.degree,...
    'diff',options.derivatives(2),'kron',true); %#ok<ASGLU>
W1=sgn*W1;
W=W1.'*wmat*W2;
y=p1.profile(:).'*W*p2.profile(:);
end
