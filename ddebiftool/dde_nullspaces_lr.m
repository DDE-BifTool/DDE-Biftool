function [v,w]=dde_nullspaces_lr(A,varargin)
%% determine left and right null space of large sparse matrix (possibly non-square)
% This routine relies on the property of the LU factorization that it puts
% the zero elements at the bottom of U
%
% $Id: dde_nullspaces_lr.m 319 2019-01-31 02:14:56Z jansieber $
%% 
default={'nulltol',1e-8,'bordered',true,'check',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% obtain first (unsafe?) approximations using LU decomposition
% for left nullspace
v=dde_nullspace(A,'check',false,'nulltol',options.nulltol,pass_on{:});
[s1,s2]=size(A);
leftnulldim=size(v,2);
%% for right nullspace
rightnulldim=s1-s2+leftnulldim;
assert(rightnulldim>=0);
w=dde_nullspace(A',pass_on{:},'nulldim',rightnulldim);
assert(rightnulldim==size(w,2));
if ~options.bordered
    return
end
%% obtain safer(?) nullspaces using bordered matrix (unnessecary?)
Aext=[A,w;v',sparse(leftnulldim,rightnulldim)];
assert(size(Aext,1)==size(Aext,2))
lrhs=cat(1,zeros(s1,leftnulldim),eye(leftnulldim));
v=Aext\lrhs;
[v,~]=qr(v(1:s2,:),0);
rrhs=cat(1,zeros(s2,rightnulldim),eye(rightnulldim));
w=Aext'\rrhs;
[w,~]=qr(w(1:s1,:),0);
if options.check
    resv=max(abs(A*v),[],1);
    resw=max(abs(A'*w),[],1);
    An=norm(A,'inf');
    errvlarge=find(resv(:)>options.nulltol*An);
    if ~isempty(errvlarge)
        [dum,ix]=max(abs(errvlarge)); %#ok<ASGLU>
        warning('dde_nullspaces_lr:residual',...
            'dde_nullspaces_lr: residual for right nullvector %d =%g\n',...
            errvlarge(ix),resv(errvlarge(ix)));
    end
    errwlarge=find(resw(:)>options.nulltol*An);
    if ~isempty(errwlarge)
        [dum,ix]=max(abs(errwlarge)); %#ok<ASGLU>
        warning('dde_nullspaces_lr:residual',...
            'dde_nullspaces_lr: residual for left nullvector %d =%g\n',...
            errwlarge(ix),resw(errwlarge(ix)));
    end    
end
end
