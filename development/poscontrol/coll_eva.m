%% Evaluate (derivative of) piecewise collocation polynomial
% assuming all requested points are inside the collocation interval [0,1]
%
%%
function [y,J]=coll_eva(profile,mesh,x,degree,varargin)
%% INPUT:
%
% * profile: profile on mesh t
% * mesh: representation points in [0,1]^(m*l+1)
% * x: point(s) where to evaluate
% * degree: order polynomials
%% OUTPUT: 
%	y value of profile at x
%   J: sparse Jacobian:  y(:)=J*profile(:)
%
%% options
default={'kron',false,'diff',0};
options=dde_set_options(default,varargin,'pass_on');
%% coarse mesh
tcoarse=mesh(1:degree:end);
nt=length(mesh);
nx=length(x);
%% find for each c the corresponding interpolation interval
% for c(i) the interval number is it(i) where it(i) points into tcoarse
[dum,itx]=sort([mesh,x]); %#ok<ASGLU>
it=NaN(1,nx);
jt=0;
jx=1;
for i=1:length(itx)
    if itx(i)<=nt
        jt=jt+1;
    else
        it(jx)=jt;
        jx=jx+1;
    end
    if jx>nx
        break
    end
end
%% undo sorting
icsort=itx(itx>nt)-nt;
it(icsort)=it;
it(x==1)=nt;
itcoarse=floor((it-1)/degree)+1;
itcoarse(itcoarse==length(tcoarse))=length(tcoarse)-1;
%% evaluate Lagrange polynomials on all interpolation times
iseq= x==mesh(it);
intlen=tcoarse(itcoarse+1)-tcoarse(itcoarse);
cscal=2*(x-tcoarse(itcoarse))./intlen-1;
%% evaluate barycentric weights
% for flexibility, (not assuming any particular interpolation grid)
[w,base_v,Dmat]=coll_barywt(mesh,degree,options.diff);
%% calculate profiles
y=NaN(size(profile,1),length(x));
ti_m=(itcoarse-1)*degree+1;
jac_ind=NaN((degree+1),nx,2);
jac_vals=zeros(degree+1,nx);
for i=1:nx
    ix=ti_m(i)+(0:degree)';
    jac_ind(:,i,1)=i;
    jac_ind(:,i,2)=ix;
    if iseq(i)
        isel=it(i)-ti_m(i)+1;
        jac_vals(isel,i)=1;
    else
        fac=w./(cscal(i)-base_v);
        jac_vals(:,i)=fac/sum(fac);
    end
    jac_vals(:,i)=Dmat'*jac_vals(:,i)/intlen(i)^options.diff;
    y(:,i)=profile(:,ix)*jac_vals(:,i);
end
%% return Jacobian if requested
if nargout>1
    jac_ind=reshape(jac_ind,[],2);
    jac_vals=jac_vals(:);
    iremove=jac_vals==0;
    jac_vals=jac_vals(~iremove);
    jac_ind=jac_ind(~iremove,:);
    J=sparse(jac_ind(:,1),jac_ind(:,2),jac_vals,nx,nt);
    if options.kron
       J=kron(J,speye(size(y,1)));
    end
end
end