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
% $Id: br_bisection.m 362 2019-07-14 15:49:40Z jansieber $
%
%%
default={'print',0,'distance',branch.method.bifurcation.secant_tolerance,...
    'min_iterations',branch.method.bifurcation.secant_iterations,...
    'stabilityfield','l0',...
    'bif_mon_reference',false};
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
    bif_monitor=@(p)std_bifs(funcs,branch.method,bif_mon_inp,p,options.stabilityfield);
else
    bif_monitor=bif_mon_inp;
end
bif_mon=@(p,pref)bif_mon_wrap(bif_monitor,p,pref,options.bif_mon_reference);
[res(1),p1]=bif_mon(p1,p1);
[res(2),p2]=bif_mon(p2,p1);
if prod(res)>0
    if options.print>0
        fprintf(['br_bisection: detection function does not change sign:',...
            'r(1)=%g, r(2)=%g'],res(1),res(2));
    end
    indbif=NaN;
    indmap=1:length(branch.point);
    notcorrected=[];
    return;
end
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
    [resnew,pnewcor]=bif_mon(pnewcor,plist(end));
    plist=[plist,pnewcor]; %#ok<AGROW>
    slist=[slist,snew]; %#ok<AGROW>
    iscorrected=[iscorrected,logical(suc)]; %#ok<AGROW>
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
plist=arrayfun(@(x)dde_trim_point(x,branch.point(inds(1))),plist);
branch.point=[branch.point(1:inds(1)),plist(2:end-1),branch.point(inds(2):end)];
[idum,indbif]=max(ix); %#ok<ASGLU>
indbif=indbif+inds(1)-1;
notcorrected=notcorrected+inds(1)-1;
end
function res=std_bifs(funcs,mth,type,p,stabilityfield)
args={'critical',true};
sargs={'args',{args}};
xargs=[{'exclude_trivial',true},args];
sxargs={'args',{xargs}};
bifs=struct(...
    'hopf', struct('res',@(lambda)real(lambda),sargs{:}),...
    'fold', struct('res',@(lambda)real(lambda),sargs{:}),...
    'pfold',struct('res',@(mu)real(mu-1),sargs{:}),...
    'tr',   struct('res',@(mu)abs(mu)-1,sxargs{:}),...
    'pd',   struct('res',@(mu)real(mu+1),sxargs{:}),...
    'hoho', struct('res',@(lambda)real(lambda),sxargs{:}),...
    'zeho', struct('res',@(lambda)real(lambda),sxargs{:}),...
    'hoze', struct('res',@(lambda)real(lambda),sxargs{:}));
[nunst,dom]=GetStability(p,'method',mth,bifs.(type).args{:},...
    'stabilityfield',stabilityfield,'funcs',funcs); %#ok<ASGLU>
res=bifs.(type).res(dom);
if isnan(res)
    res=-Inf;
end
end
%% wrapper around bif_mon
function [res,pnew]=bif_mon_wrap(bif_mon,p,pref,hasreference)
if hasreference
    [res,pnew]=bif_mon(p,pref);
else
    pnew=p;
    res=bif_mon(p);
end
end
