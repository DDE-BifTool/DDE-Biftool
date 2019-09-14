function [mu,eigenfuncs]=mult_app(funcs,point,rho,max_number,col,d_ac)
%% find Floquet multipliers & modes of periodic orbit
% function [mu,eigenfuncs]=mult_app(period,profile,mesh,degree,rho,max_number,col,par,d_ac)
% INPUT: 
%   funcs problem functions
%   point of solution type psol
%	max_number keep at most max_number multipliers 
%	col collocation points
%       par current parameter values in R^p
%       d_ac (only for state-dependent delays) tau<d_ac is treated as 
%             tau<0 (stability is not computed)
% OUTPUT:
%       mu approximations of requested multipliers
%       eigenfuncs (if requested) corresponding modes

%  (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: mult_app.m 27 2014-04-14 18:38:00Z jan.sieber $
%
%%
profile=point.profile;
if nargin<10
    d_ac=1e-8;
end
%% if mesh is empty assume equidistant mesh:
%% obtain Jacobian
[J,resdum,tT,extmesh]=psol_jac_reduced(funcs,point,[],'wrapJ',false); %#ok<ASGLU>
n=size(profile,1);
delays=tT(2:end,:);
if nargout>1
    geteigenfuncs=true;
else
    geteigenfuncs=false;
end    
if numel(delays)>0 && min(delays(:))<d_ac
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
    M0=-J(:,n_ext+1:end)\J(:,1:n_ext);
    if n_ext<=s1
        M=M0(end-n_ext+1:end,:);
    else
        M=[zeros(n_ext-s1,s1),eye(n_ext-s1);M0];
    end
    if isempty(M)
        mu=[];
        return;
    end
    M=full(M);
    if ~geteigenfuncs
        s=eig(M);
    else
        [ef,s]=eig(M);
        s=diag(s);
    end
else
    %% negative delays present or eigenvectors required
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
