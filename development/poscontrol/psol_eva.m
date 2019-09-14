%% Evaluate piecewise periodic collocation polynomial
function [px,J]=psol_eva(pt,x,varargin)
%% function [px,J]=psol_eva(profile1,t,x,m);
% INPUT:
% * pt with profile: profile on mesh t, mesh: representation points in
% [0,1]^(m*l+1), degree: order polynomials
% * x: point(s) where to evaluate
% OUTPUT: 
%	px value of profile at x
%   J: sparse Jacobian:  y(:)=J*profile(:) if wrapJ is true, otherwise
%   J is structure with row indices, col indices and values for sparse
%   matrix, but which may contain negative column indices
%  
%
% (c) DDE-BIFTOOL v. 3.1.1(106), 22/08/2015
%
%%
default={'wrapJ',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% coarse mesh
%% wrap requested x into [0,1)
c=x-floor(x);
c(x==1)=1; % for rotations
[px,J]=coll_eva(pt.profile,pt.mesh,c,pt.degree,'kron',true,pass_on{:});
if ~options.wrapJ
    %% shift column indices toward past if x was negative
    n=size(pt.profile,1);
    nt=length(pt.mesh)-1;
    ishift=floor(x(:));
    sel1=x==floor(x)&x>0;
    ishift(sel1)=ishift(sel1)-1;
    ishiftn=reshape(repmat(ishift,1,n)',[],1);
    [ir,ic,Jvals]=find(J);
    %% determine column indices of points outside base interval
    ic=ic+ishiftn(ir)*n*nt;
    imin=min(ishift)*n*nt+1; % minimal index (whole base intervals)
    subint=floor((min(ic)-imin)/pt.degree);
    imin=imin+subint*pt.degree;
    %% store structure information for Jacobian
    J=struct('ir',ir,'ic',ic,'vals',Jvals,...
        'nr',size(J,1),'nc',size(J,2),'imin',imin,'x',x,'nkron',n);
end
end