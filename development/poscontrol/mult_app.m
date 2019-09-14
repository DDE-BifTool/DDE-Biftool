function [mu,eigenfuncs]=mult_app(funcs,point,rho,max_number,varargin)
%% find Floquet multipliers & modes of periodic orbit
% function [mu,eigenfuncs]=mult_app(period,profile,mesh,degree,rho,max_number,col,par)
% INPUT: 
%   funcs problem functions
%	point with profile, mesh, parameter, degree
%   rho keep multipliers with modulus >= rho
%	max_number keep at most max_number multipliers 
%   (named option) d_ac (only for state-dependent delays) tau<d_ac is treated as 
%             tau<0 (stability is not computed)
%   (named option) sparse (default false) use eigs if true
%   all other parameters get passed on to psol_jac_sparse
% OUTPUT:
%       mu approximations of requested multipliers
%       eigenfuncs (if requested) corresponding modes

%  (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: mult_app.m 27 2014-04-14 18:38:00Z jan.sieber $
%
%%
default={'d_ac',1e-8,'sparse',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% if mesh is empty assume equidistant mesh:
if isempty(point.mesh)
  point.mesh=0:1/(size(profile,2)-1):1;
end;
%% obtain Jacobian
if options.sparse
    [J,resdum,tT,extmesh]=psol_jac_combined(funcs,point,[],pass_on{:},...
        'wrapJ',false,'period',false,'bc',false,'ph',false); %#ok<ASGLU>
    n=size(point.profile,1);
    delays=tT(:,2:end);
else
    [J,resdum,tT,extmesh]=psol_jac(funcs,[],point.period,point.profile,...
        point.mesh,point.degree,point.parameter,...
        [],false,'wrapJ',false); %#ok<ASGLU>
    n=size(point.profile,1);
    delays=tT(2:end,:);
end
if nargout>1
    geteigenfuncs=true;
else
    geteigenfuncs=false;
end    
if numel(delays)>0 && min(delays(:))<options.d_ac
    % solve full eigenvalue problem
    solvefull=true;
else
    solvefull=false;
end
if ~solvefull
    %% DDE or ODE compute monodromy matrix
    % was: M=monodromy_matrix(J,n,m,extmesh);
    [s1,s2]=size(J);
    n_ext=s2-s1;
    Mrhs=-J(:,1:n_ext);
    if ~options.sparse
        M0=J(:,n_ext+1:end)\Mrhs;
        if n_ext<=s1
            M=M0(end-n_ext+1:end,:);
        else
            M=[zeros(n_ext-s1,s1),eye(n_ext-s1);M0];
        end
        if isempty(M)
            mu=[];
            return
        end
        if ~geteigenfuncs
            s=eig(full(M));
        else
            [ef,s]=eig(full(M));
            s=diag(s);
        end
    else %% use eigs
        [L,U,P,Q,R]=lu(J(:,n_ext+1:end));
        Mfun=@(x)Monodromy_func(x,Q,U,L,P,R,Mrhs);
        sMhrs=size(Mrhs,2);
        if ~geteigenfuncs
            s=eigs(Mfun,sMhrs,min(max_number,sMhrs-1));
        else
            [ef,s]=eigs(Mfun,sMhrs,min(max_number,sMhrs-1));
            s=diag(s);
        end
        
    end
else % no sparse method here
    %% negative delays present
    ll=length(extmesh);
    n_ext=sum(extmesh<=0|extmesh>1)*n;
    B=zeros(ll*n);
    B(end-n_ext+1:end,1:n_ext)=eye(n_ext);
    J=[full(J);...
        zeros(n_ext,ll*n-n_ext),eye(n_ext)];
    if ~geteigenfuncs
        s=eig(J,B);
    else
        [ef,s]=eig(J,B);
        s=diag(s);
    end
    sel=~isinf(s) & ~isnan(s);
    s=s(sel);
    if geteigenfuncs
        ef=ef(:,sel);
    end
end
[dummy,I]=sort(abs(s)); %#ok<ASGLU>

mu=s(I(end:-1:1));
mu=mu(1:min(max_number,length(mu)));
sel=abs(mu)>=rho;
mu=mu(sel);
if geteigenfuncs
    ef=ef(:,I(end:-1:1));
    ef=ef(:,1:min(max_number,length(mu)));
    ef=ef(:,sel);
    if ~solvefull
        n_ext=size(J,2)-size(J,1);
        dim=size(profile,1);
        ef1=-J(:,n_ext+1:end)\J(:,1:n_ext)*ef;
        eigenfuncs=[ef(end-dim+1:end,:);ef1];
    else
        eigenfuncs=ef;
    end
end
end
%% function for sparse eigs
function Mx=Monodromy_func(x,Q,U,L,P,R,Mrhs)
Mx1=Q*(U\(L\(P*(R\(Mrhs*x)))));
if size(Mx1,1)<size(x,1)
    Mx=[x(size(Mx1,1)+1:end,:);Mx1];
else
    Mx=Mx1(size(Mx1,1)-size(x,1)+1:end,:);
end
end