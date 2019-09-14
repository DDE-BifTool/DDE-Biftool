function J=dde_gen_deriv(dffunc,xx,par,nx,np,v)
%% compose derivatives of r.h.s wrt state and parameter from directional derivatives
%
% function J=dde_gen_f_deriv(dirderi,xx,par,nx,np,v)
%
%% INPUT:
%   dffunc: directional derivative function df(order,x,p,dx,dp) 
%	xx state variable and delayed state variables columnwise
%	par list of parameter values
%	nx empty or list of requested state-derivatives (numbers of delay or zero)
%	np empty or list of requested parameter-derivatives
%	v matrix to multiply result with for 2nd xx derivative
%% OUTPUT:
%	J result of derivatives on righthandside multiplied with v
%% COMMENT:
%	pass on numerical derivatives for dffunc if necessary
%
% <html>
% $Id: dde_gen_deriv.m 174 2017-03-10 21:59:17Z jansieber $
% </html>
%
%%
xsize=size(xx);
n=xsize(1);
ndelays=xsize(2);
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xx=reshape(xx,[n,ndelays,nvec]);

if length(nx)==1 && isempty(np) && isempty(v)
    %% first order derivatives of the state:
    xdev=zeros(n,ndelays,n);
    xdev(:,nx+1,:)=reshape(eye(n),[n,1,n]);
    xdev=repmat(xdev,[1,1,1,nvec]);
    xx=repmat(reshape(xx,[n,ndelays,1,nvec]),[1,1,n,1]);
    xx=reshape(xx,[n,ndelays,n*nvec]);
    xdev=reshape(xdev,size(xx));
    df=dffunc(1,xx,par,xdev,zeros(size(par)));
    J=reshape(df,[size(df,1),n,vecdim]);
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% first order derivative with respect to parameter
    dpar=zeros(size(par));
    dpar(np)=1;
    J=dffunc(1,xx,par,zeros(size(xx)),dpar);
    J=reshape(J,size(J,1),1,nvec);
elseif length(nx)==2 && isempty(np) && ~isempty(v)
    %% second order state derivatives, applied to deviation
    xx=repmat(reshape(xx,[n,ndelays,1,nvec]),[1,1,n,1]);
    xx=reshape(xx,[n,ndelays,n*nvec]);
    dev0={v(:,ones(1,n)),eye(n)};
    for i=2:-1:1
        xdev{i}=zeros(n,ndelays,n);
        xdev{i}(:,nx(i)+1,:)=reshape(dev0{i},[n,1,n]);
        xdev{i}=repmat(xdev{i},[1,1,1,nvec]);
        xdev{i}=reshape(xdev{i},size(xx));
    end
    df{1}=dffunc(2,xx,par,xdev{1}+xdev{2},zeros(size(par)));
    df{2}=dffunc(2,xx,par,xdev{1}-xdev{2},zeros(size(par)));
    J=0.25*(df{1}-df{2});
    J=reshape(J,[size(J,1),n,vecdim]);
elseif length(nx)==2 && isempty(np) && isempty(v)
    %% complete second order state derivatives (relevant for delays)
    % not sure this is used anywhere currently
    xx=repmat(reshape(xx,[n,ndelays,1,1,nvec]),[1,1,n,n,1]);
    xx=reshape(xx,[n,ndelays,n*n*nvec]);
    e1=repmat(eye(n),[1,1,n]);
    e2=repmat(reshape(eye(n),[n,1,n]),[1,n,1]);
    dev0={e1,e2};
    for i=2:-1:1
        xdev{i}=zeros(n,ndelays,n,n);
        xdev{i}(:,nx(i)+1,:,:)=reshape(dev0{i},[n,1,n,n]);
        xdev{i}=repmat(xdev{i},[1,1,1,1,nvec]);
        xdev{i}=reshape(xdev{i},size(xx));
    end
    df{1}=dffunc(2,xx,par,xdev{1}+xdev{2},zeros(size(par)));
    df{2}=dffunc(2,xx,par,xdev{1}-xdev{2},zeros(size(par)));
    J=0.25*(df{1}-df{2});
    assert(size(J,1)==1);
    J=reshape(J,[n,n,vecdim]);
elseif length(nx)==1 && length(np)==1 && isempty(v)
    %% mixed state-parameter derivatives:
    xdev=zeros(n,ndelays,n);
    xdev(:,nx+1,:)=reshape(eye(n),[n,1,n]);
    xdev=repmat(xdev,[1,1,1,nvec]);
    dpar=zeros(size(par));
    dpar(np)=1;
    xx=repmat(reshape(xx,[n,ndelays,1,nvec]),[1,1,n,1]);
    xx=reshape(xx,[n,ndelays,n*nvec]);
    xdev=reshape(xdev,size(xx));
    df{1}=dffunc(2,xx,par,xdev,dpar);
    df{2}=dffunc(2,xx,par,xdev,-dpar);
    J=0.25*(df{1}-df{2});
    J=reshape(J,[size(J,1),n,vecdim]);
elseif length(nx)==1 && isempty(np) && ~isempty(v)
    %% higher-order state derivatives in direction v
    xdev=reshape(v,size(xx,1),size(xx,2),[]);
    df=dffunc(nx,xx,par,xdev,zeros(size(par)));
    J=reshape(df,[size(df,1),n,vecdim]);
else
    %% not implemented
    sv=num2str(size(v));
    error('dde_gen_deriv:notimplented',...
        ['dde_gen_deriv: requested derivative nx=%d, np=%d, size(v)=(%s)',...
        ' does not exist!'],nx,np,sv);
end
end
