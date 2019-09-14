function J=dde_gen_deriv(dffunc,xx,par,nx,np,v,directional_derivative)
%% compose derivatives of r.h.s wrt state and parameter from directional derivatives
%
% function J=dde_gen_deriv(dirderi,xx,par,nx,np,v)
%
%% INPUT:
%   dffunc: directional derivative function df(order,x,p,dx,dp)
%	xx state variable and delayed state variables columnwise
%	par list of parameter values
%	nx empty or list of requested state-derivatives (numbers of delay or zero)
%	np empty or list of requested parameter-derivatives
%	v matrix to multiply result with for 2nd xx derivative
%   directional_derivative boolean whether directional derivatives or
%   multilinear form files are used
%% OUTPUT:
%	J result of derivatives on righthandside multiplied with v
%% COMMENT:
%	pass on numerical derivatives for dffunc if necessary
%
% <html>
% $Id: dde_gen_deriv.m 348 2019-06-19 13:09:01Z jansieber $
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
    xdevr=reshape(xdev,size(xx));
    if size(par,3)==nvec
        par=repmat(par,1,1,n);
    end
    df=dffunc{1}(xx,par,xdevr,zeros(size(par)));
    J=reshape(df,[size(df,1),n,vecdim]);
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% first order derivative with respect to parameter
    dpar1=zeros(size(par));
    dpar1(1,np,:)=1;
    J=dffunc{1}(xx,par,zeros(size(xx)),dpar1);
    J=reshape(J,[size(J,1),1,vecdim]);
elseif length(nx)==2 && isempty(np) && ~isempty(v)
    %% second order state derivatives, applied to deviation
    if directional_derivative
        xx=repmat(reshape(xx,[n,ndelays,1,nvec]),[1,1,n,1]);
        xx=reshape(xx,[n,ndelays,n*nvec]);
        v=reshape(v,n,1,nvec);
        dev0={v(:,ones(1,n),:),repmat(eye(n),[1,1,nvec])};
        for i=2:-1:1
            xdev{i}=zeros(n,ndelays,n,nvec);
            xdev{i}(:,nx(i)+1,:,:)=reshape(dev0{i},[n,1,n,nvec]);
            xdev{i}=reshape(xdev{i},size(xx));
        end
        if size(par,3)==nvec
            par=repmat(par,1,1,n);
        end
        df{1}=dffunc{2}(xx,par,xdev{1}+xdev{2},zeros(size(par)));
        df{2}=dffunc{2}(xx,par,xdev{1}-xdev{2},zeros(size(par)));
        J=0.25*(df{1}-df{2});
    else
        xx=repmat(reshape(xx,[n,ndelays,1,nvec]),[1,1,n,1]);
        xx=reshape(xx,[n,ndelays,n*nvec]);
        xdev=zeros(n,ndelays,2,n,nvec);
        v=reshape(v,n,1,1,1,nvec);
        xdev(:,nx(1)+1,1,:,:)=v(:,:,:,ones(1,n),:);
        xdev(:,nx(2)+1,2,:,:)=...
            reshape(repmat(eye(n),[1,nvec]),[n,1,n,1,1,nvec]);
        np=length(par);
        pdev=zeros(np,2,n,nvec);
        par=repmat(par,1,1,n,nvec);
        J=dffunc{2}(xx,par,xdev,pdev);
    end
    J=reshape(J,[size(J,1),n,vecdim]);
    % elseif length(nx)==2 && isempty(np) && isempty(v)
    %     %% complete second order state derivatives (relevant for delays)
    %     % not sure this is used anywhere currently
    %     xx=repmat(reshape(xx,[n,ndelays,1,1,nvec]),[1,1,n,n,1]);
    %     xx=reshape(xx,[n,ndelays,n*n*nvec]);
    %     e1=repmat(eye(n),[1,1,n]);
    %     e2=repmat(reshape(eye(n),[n,1,n]),[1,n,1]);
    %     dev0={e1,e2};
    %     for i=2:-1:1
    %         xdev{i}=zeros(n,ndelays,n,n);
    %         xdev{i}(:,nx(i)+1,:,:)=reshape(dev0{i},[n,1,n,n]);
    %         xdev{i}=repmat(xdev{i},[1,1,1,1,nvec]);
    %         xdev{i}=reshape(xdev{i},size(xx));
    %     end
    %     df{1}=dffunc{2}(xx,par,xdev{1}+xdev{2},zeros(size(par)));
    %     df{2}=dffunc{2}(xx,par,xdev{1}-xdev{2},zeros(size(par)));
    %     J=0.25*(df{1}-df{2});
    %     assert(size(J,1)==1);
    %     J=reshape(J,[n,n,vecdim]);
elseif length(nx)==1 && length(np)==1 && isempty(v)
    %% mixed state-parameter derivatives:
    if directional_derivative
        xdev=zeros(n,ndelays,n);
        xdev(:,nx+1,:)=reshape(eye(n),[n,1,n]);
        xdev=repmat(xdev,[1,1,1,nvec]);
        if size(par,3)==nvec
            par=repmat(par,1,1,n);
        end
        dpar1=zeros(size(par));
        dpar1(1,np,:)=1;
        xx=repmat(reshape(xx,[n,ndelays,1,nvec]),[1,1,n,1]);
        xx=reshape(xx,[n,ndelays,n*nvec]);
        xdev=reshape(xdev,size(xx));
        df{1}=dffunc{2}(xx,par,xdev,dpar1);
        df{2}=dffunc{2}(xx,par,xdev,-dpar1);
        J=0.25*(df{1}-df{2});
    else
        xx=repmat(reshape(xx,[n,ndelays,1,nvec]),[1,1,n,1]);
        xx=reshape(xx,[n,ndelays,n*nvec]);
        xdev=zeros(n,ndelays,2,n,nvec);
        xdev(:,nx+1,1,:,:)=...
            reshape(repmat(eye(n),[1,nvec]),[n,1,n,1,1,nvec]);
        pdev=zeros(length(par),2,n,nvec);
        pdev(np,2,:,:)=1;
        par=repmat(par,1,1,n,nvec);
        J=dffunc{2}(xx,par,xdev,pdev);       
    end
    J=reshape(J,[size(J,1),n,vecdim]);
elseif length(nx)==1 && isempty(np) && ~isempty(v)
    %% higher-order state derivatives in direction v
    xdev=reshape(v,size(xx,1),size(xx,2),[]);
    df=dffunc{nx}(xx,par,xdev,zeros(size(par)));
    J=reshape(df,[size(df,1),n,vecdim]);
elseif isempty(nx) && length(np)==2
    %% 2nd-order derivative with respect to parameter
    dpar1=zeros(size(par));
    dpar1(1,np(1),:)=1;
    dpar2=zeros(size(par));
    dpar2(1,np(2),:)=1;
    J1=dffunc{2}(xx,par,zeros(size(xx)),dpar1+dpar2);
    J2=dffunc{2}(xx,par,zeros(size(xx)),dpar1-dpar2);
    J=0.25*(J1-J2);
    J=reshape(J,[size(J,1),1,vecdim]);
else
    %% not implemented
    sv=num2str(size(v));
    error('dde_gen_deriv:notimplented',...
        ['dde_gen_deriv: requested derivative nx=%d, np=%d, size(v)=(%s)',...
        ' does not exist!'],nx,np,sv);
end
end
