function [px,it,cscal]=psol_eva(profile1,t,x,m)

% function [px,it,cscal]=psol_eva(profile1,t,x,m);
% INPUT:
%	profile1 profile on mesh t
%	t representation points in [0,1]^(m*l+1)
%	x point(s) where to evaluate
%	m order polynomials
% OUTPUT: 
%	px value of profile at x
%   it mesh intervals where x was located
%   cscal interpolation point inside mesh (scaled up by interval length)
%
% (c) DDE-BIFTOOL v. 1.00, 08/04/2000
%
% $Id: psol_eva.m 106 2015-08-22 20:25:33Z jansieber $
%
%% coarse mesh
tcoarse=t(1:m:end);
nt=length(tcoarse);
nx=length(x);
%% wrap requested x into [0,1)
c=x-floor(x);
%% find for each c the corresponding interpolation interval
% for c(i) the interval number is it(i) where it(i) points into tcoarse
[tdum,itx]=sort([tcoarse,c]); %#ok<ASGLU>
it=NaN(1,nx);
jt=0;
jx=1;
for i=1:length(itx)
    if itx(i)<=nt
        jt=jt+1;
    else
        it(jx)=jt;
        jx=jx+1;
    end
    if jx>nx
        break
    end
end
%% undo sorting
icsort=itx(itx>nt)-nt;
it(icsort)=it;
it(x==1)=nt-1;
c(x==1)=1;
%% evaluate Lagrange polynomials on all interpolation times
cscal=(c-tcoarse(it))./(tcoarse(it+1)-tcoarse(it));
Pa=poly_elg(m,cscal);
if length(x)==1
    Pa=Pa(:);
end
%% calculate profiles
px=NaN(size(profile1,1),length(x));
ti_m=(it-1)*m+1;
for i=1:length(x)
    ix=ti_m(i)+(0:m);
    px(:,i)=profile1(:,ix)*Pa(:,i);
end
end
%% old version (shorter but slower for p_remesh etc)
% for i=1:length(x)
%     
%     c=x(i)-floor(x(i));
%     
%     index=length(t)-m;
%     while (c<t(index))
%         index=index-m;
%     end;
%     Pa=poly_elg(m,(c-t(index))/(t(index+m)-t(index)));
%     px(:,i)=profile1(:,index:index+m)*Pa';
% end
