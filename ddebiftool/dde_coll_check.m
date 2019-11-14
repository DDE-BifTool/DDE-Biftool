function [pt,submesh]=dde_coll_check(pt)
%% Check if profile, degree and mesh are compatible, and insert mesh
%
% $Id: dde_coll_check.m 369 2019-08-27 00:07:02Z jansieber $
%%
assert(~isempty(pt.degree),'dde_coll_check:mesh','dde_coll_check: polynomial degree needs to be specified');
ntst=size(pt.profile,2);
assert(mod((ntst-1)/pt.degree,1)==0,...
    'dde_coll_check:mesh',...
    'dde_coll_check: degree %d and size of profile (%d,%d) incompatible',...
    pt.degree,size(pt.profile,1),size(pt.profile,2));
if isempty(pt.mesh)
    pt.mesh=linspace(0,1,ntst);
else
    assert(length(pt.mesh)==size(pt.profile,2),...
        'dde_coll_check:mesh',...
        'dde_coll_check: mesh length %d and size of profile (%d,%d) incompatible',...
        length(pt.mesh),size(pt.profile,1),size(pt.profile,2));
    tcoarse=pt.mesh(1:pt.degree:end);
    hcoarse=diff(tcoarse);
    od=ones(pt.degree,1);
    hmesh=reshape(pt.mesh(2:end),pt.degree,[]);
    hmesh=(hmesh-tcoarse(od,1:end-1))./hcoarse(od,:);
    hmeshd=diff(hmesh,[],2);
    assert(isempty(hmeshd)||max(abs(hmeshd(:)))<1e-8,'dde_coll_check:mesh',...
        'dde_coll_check: internal meshes differ');
end
submesh=pt.mesh(1:pt.degree+1)/pt.mesh(pt.degree+1);
end
