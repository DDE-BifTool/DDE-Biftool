function v=dde_nullspace(A,varargin)
%% determine kernel of large sparse matrix (possibly non-square)
% This routine relies on the property of the LU factorization that it puts
% the zero elements at the bottom of U
%
% $Id: dde_nullspace.m 319 2019-01-31 02:14:56Z jansieber $
%% 
default={'nulltol',1e-8,'check',true,'nulldim',[]};
options=dde_set_options(default,varargin,'pass_on');
[s1,s2]=size(A);
if s1<s2
    A=cat(1,A,sparse(s2-s1,s2));
end
[L,U,P,Q,D]=lu(sparse(A)'); %#ok<ASGLU>
%% determine presumable kernel dimension
if isempty(options.nulldim)
    uind=find(full(max(abs(U),[],2).'<options.nulltol));
else
    [dum,iu]=sort(full(max(abs(U),[],2)),'descend'); %#ok<ASGLU>
    uind=iu(end-options.nulldim+1:end);
end
dim=length(uind);
w=full(sparse(uind(:),(1:dim),ones(dim,1),s2,dim));
v=P'*((L')\w);
d=diag(D);
v=v./d(:,ones(size(v,2),1));
[v,~]=qr(full(v),0);
if options.check
    res=max(abs(A*v),[],1);
    An=norm(A,'inf');
    errlarge=find(res(:)>options.nulltol*An);
    if ~isempty(errlarge)
        [dum,ix]=max(abs(errlarge)); %#ok<ASGLU>
        warning('dde_nullspace:residual',...
            'dde_nullspace: residual for nullvector %d =%g\n',errlarge(ix),res(errlarge(ix)));
    end
end
end
