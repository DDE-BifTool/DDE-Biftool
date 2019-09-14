function point=p_remesh(point,inp2,new_mesh,varargin)
%% Remesh collocation mesh to new given mesh or to equidistribute error
%
% function rm_point=p_remesh(point,new_degree,new_mesh)
% INPUT:
%	point: solution point (with fields mesh and profile)
%	(optional) inp2: new degree for new representation or poitn from which
%	mesh is copied
%	(optional) new_mesh: new mesh or new number of intervals
% OUTPUT:
%	rm_point interpolated point on new mesh
%
% $Id: p_remesh.m 369 2019-08-27 00:07:02Z jansieber $
%
% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%% routine is trivial for points without mesh or with single interval
if ~isfield(point,'mesh') || ~isfield(point,'degree') ||...
        length(point.mesh)<=point.degree+1
    return
end
if nargin<2
    %% with single argument reset mesh to equidistribute error
    new_degree=point.degree;
    new_mesh=length(point.mesh(1:point.degree:end));
elseif nargin<3
    %% with two arguments: 
    %% if second argument is also point, copy mesh
    new_mesh=length(point.mesh(1:point.degree:end));
    if isstruct(inp2) && isfield(inp2,'degree') && isfield(inp2,'mesh')
        inp2=dde_coll_check(inp2);
        new_mesh=inp2.mesh;
        new_degree=inp2.degree;
    else
        %% otherwise, 2nd argument is new degree
        new_degree=inp2;
    end
else
    new_degree=inp2;
end
if length(new_mesh)==length(point.mesh) && all(point.mesh==new_mesh) &&...
        new_degree==point.degree
    return
end
[point,conversion]=dde_coll_convert(point);
[point,submesh]=dde_coll_check(point);
if length(new_mesh)==1
    new_mesh=dde_coll_newmesh(point,new_mesh,new_degree,'grid',submesh,varargin{:});
end
point.profile=dde_coll_eva(point.profile,point.mesh,new_mesh,point.degree,varargin{:});
point.mesh=new_mesh;
point.degree=new_degree;
point=dde_coll_check(point);
point=dde_coll_convert(point,conversion);
if isfield(point,'stability')
    point.stability=[];
end
end
