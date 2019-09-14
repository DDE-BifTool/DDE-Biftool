function [Marg,MPfun]=dde_psol_monodromy(J,closest,matrix)
%% return monodromy matrix in functional form from collocation Jacobian
%
% Input
%
% * J: jacobian from dde_coll_jac extended into the past as far as
% necessary
% * closest (optional): if present create function for (M-c)^(-1)
%
% Output
%
% * Marg: cell of first one (or two) args of eig(s). Cases: sparse & lm:
% {Mfun,n_ext} where Mfun applies monodromy matrix to vector x of length
% n_ext, returns vector of length n_ext. full: monodromy matrix {A}.
% sparse & closest: sparse matrix pair {A,B}
% * n_ext: length of history vector
% * MPfun(x): applies monodromy matrix to vector x of length n_ext, returns
% (longer) vector of length size(J,2)
%
% $Id: dde_psol_monodromy.m 350 2019-06-19 13:22:44Z jansieber $
%%
if nargin<3 || strcmp(matrix,'full')
    M_is_sparse=false;
elseif strcmp(matrix,'sparse')
    M_is_sparse=true;
end
[s1,s2]=size(J);
n_ext=s2-s1;
if nargin>1 && ~isempty(closest) && M_is_sparse
    O=sparse(n_ext,s1);
    id=speye(n_ext);
    MA=[O,id;J];
    MB=[id,O;sparse(s1,s2)];
    Marg={MA,MB};
    MPfun=@(x)x;
else 
    n2_ext=max(n_ext-s1,0);
    B=blkdiag(speye(n2_ext),J(1:s1,n_ext+1:s2));
    Atop=cat(2,sparse(n2_ext,s1),speye(n2_ext));
    A=cat(1,Atop(:,1:n_ext),-J(1:s1,1:n_ext)); % s1 x n_ext
    Ptau=cat(2,sparse(n_ext,max(0,2*s1-s2)),speye(n_ext));
    [Lm,Um,Pm,Qm,Rm]=lu(B);
    Mfun=@(x)Ptau*(Qm*(Um\(Lm\(Pm*(Rm\(A*x))))));
    if ~M_is_sparse
        Marg={Mfun(eye(n_ext))};
    else
        Marg={Mfun,n_ext};
    end
    MPfun=@(x)cat(1,x,Qm(end-s1+1:end,:)*(Um\(Lm\(Pm*(Rm\(A*x))))));
end
end
