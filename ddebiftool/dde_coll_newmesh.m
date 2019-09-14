function t_new=dde_coll_newmesh(point,l_t_new,m_new,varargin)
%% create new mesh equidistributing the error
% function [t_new]=dde_coll_newmesh(point,l_t_new,m_new,...);
% INPUT:
%   point: solution profile
%	l_t_new size of new mesh (number of intervals)
%	m_new new number of collocation points (per interval)
% OUTPUT:
%	t_new new, adapted mesh

% (c) DDE-BIFTOOL v. 1.00, 11/03/2000
%
% $Id: dde_coll_newmesh.m 369 2019-08-27 00:07:02Z jansieber $
%
%% re-distribute coarse mesh
eqf=auto_eqd(point,'new_degree',m_new,varargin{:});

%% uniformly divide the range of eqf:
uneq=linspace(0,eqf(end),l_t_new);
ti_new=interp1(eqf,point.mesh(1:point.degree:end),uneq,'linear');
%% insert internal storage points of collocation polynomial
t_new=dde_coll_meshfill(ti_new,m_new,varargin{:});
end
