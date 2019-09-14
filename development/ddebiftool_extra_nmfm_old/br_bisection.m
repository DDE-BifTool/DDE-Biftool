function [branch,indbif,indmap,notcorrected]=...
    br_bisection(funcs,branch,inds,bif_mon_inp,varargin)
%% refine branch with bisection to improve approximate biftype
% function [refined_branch,indc,indmap]=br_bisection(funcs,branch,inds,bif)
% INPUT:
%   funcs problem functions
%	branch branch of points
%	inds indices in point list on both sides of bifurcation
%	bif_mon: bifurcation monitoring function res=bif_mon(point) for point
%	from branch, can be string 'fold', 'hopf', 'pfold','tr','pd'
% OUTPUT:
%	refined_branch branch with refinements
%   indbif index closest to bifurcation
%   indmap point with index i in original branch, is mapped to indmap(i)
%   notcorrected: list of indices in branch.point, where correction has
%   failed
%
% $Id: br_bisection.m 109 2015-08-31 23:45:11Z jansieber $
%
%%
default={'print',0,'distance',branch.method.bifurcation.secant_tolerance,...
    'min_iterations',branch.method.bifurcation.secant_iterations,'stabilityfield','l1'};
options=dde_set_options(default,varargin,'pass_on');
%% list of possible bifurcations
%% setup bisection
free_par=branch.parameter.free;
p1=branch.point(inds(1));
p2=branch.point(inds(2));
secant=p_axpy(-1,p1,p2);
secant_normalized=p_secant(secant,p_norm(p1));
psolve=@(p)p_correc(funcs,p,free_par,secant_normalized,branch.method.point,0,p1);
pstab=@(p)p_stabil(funcs,p,branch.method.stability);
ppred=@(s)p_axpy(s,secant,p1);
if ischar(bif_mon_inp)
    bif_mon=@(p)std_bifs(branch.method,bif_mon_inp,p,options.stabilityfield);
else
    bif_mon=bif_mon_inp;
end
res=[bif_mon(p1),bif_mon(p2)];
assert(prod(res)<=0,'refin:signchange',...
    'detection function does not change sign: r(1)=%g, r(2)=%g',res(1),res(2));
if options.print>0
    fprintf('Bisection residual: res(1)=%g, res(2)=%g\n',res(1),res(2));
end
dist_ini=p_norm(secant);
dist=dist_ini;
absres=min(abs(res));
s=[0,1];
plist=[p1,p2];
slist=s;
iscorrected=true(1,2);
it=0;
while (it<options.min_iterations || dist>options.distance) && absres>0
    it=it+1;
    snew=mean(s);
    pnewpred=ppred(snew);
    [pnewcor,suc]=psolve(pnewpred);
    if ~suc
        %% don't correct if correction fails
        pnewcor=pnewpred;
    end
    if isfield(plist(1),'stability')
        pnewcor.stability=pstab(pnewcor);
    end
    plist=[plist,pnewcor]; %#ok<AGROW>
    slist=[slist,snew]; %#ok<AGROW>
    iscorrected=[iscorrected,logical(suc)]; %#ok<AGROW>
    resnew=bif_mon(pnewcor);
    if options.print>0
        fprintf('Bisection: new residual=%g\n',resnew);
    end
    if resnew*res(1)<0
        s=[s(1),snew];
        res=[res(1),resnew];
    elseif resnew*res(2)<0
        s=[snew,s(2)];
        res=[resnew,res(2)];
    end
    absres=min(abs(res));
    dist=abs(diff(slist(end-1:end)))*dist_ini;
end
[sdum,ix]=sort(slist); %#ok<ASGLU>
plist=plist(ix);
notcorrected=find(~iscorrected(ix));
indmap=[1:inds(1),length(plist)-2+(inds(2):length(branch.point))];
branch.point=[branch.point(1:inds(1)),plist(2:end-1),branch.point(inds(2):end)];
[idum,indbif]=max(ix); %#ok<ASGLU>
indbif=indbif+inds(1)-1;
notcorrected=notcorrected+inds(1)-1;
end
function res=std_bifs(mth,type,p,stabilityfield)
bifs=struct(...
    'hopf',struct('res',@(lambda)real(lambda),'args',{{'critical',true}}),...
    'fold',struct('res',@(lambda)real(lambda),'args',{{'critical',true}}),...
    'pfold',struct('res',@(mu)real(mu-1),'args',{{'exclude_trivial',true,'critical',true}}),...
    'tr',struct('res',@(mu)abs(mu)-1,'args',{{'exclude_trivial',true,'critical',true}}),...
    'pd',struct('res',@(mu)real(mu+1),'args',{{'exclude_trivial',true,'critical',true}}),...
    'hoho',struct(...
     'res',@(lambda)real(lambda),...
     'args',{{...
      'critical',true,...
      'exclude_trivial',true,...
      'locate_trivial',@(p)[1i*p.omega,-1i*p.omega]}}),...
    'zeho',struct(...
     'res',@(lambda)real(lambda),...
     'args',{{...
      'critical',true,...
      'exclude_trivial',true,...
      'locate_trivial',@(p)[1i*p.omega,-1i*p.omega]}}),...
    'hoze',struct(...
     'res',@(lambda)real(lambda),...
     'args',{{...
      'critical',true,...
      'exclude_trivial',true,...
      'locate_trivial',@(p)0}}));
[nunst,dom]=GetStability(p,'method',mth,bifs.(type).args{:},...
    'stabilityfield',stabilityfield); %#ok<ASGLU>
res=bifs.(type).res(dom);
if isnan(res)
    res=-Inf;
end
end