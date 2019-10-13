function [hcli,stst]=dde_hcli_from_psol(psol,varargin)
%% convert periodic orbit to connecting orbit
default={'connections',1,'eqtime_ind',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% double psol
ntst=length(psol.mesh)-1;
psol=dde_psol_create('point',psol,...
    'mesh',[psol.mesh/2,0.5+psol.mesh(2:end)/2],...
    'period',2*psol.period,...
    'profile',[psol.profile,psol.profile(:,2:end)]);
if isempty(options.eqtime_ind)
    %% find point(s) w min derivative
    npd=sqrt(sum(dde_coll_eva(psol.profile,psol.mesh,psol.mesh,psol.degree,...
        'diff',1).^2,1));
    [dum, pos]=min(npd(1:ntst-1)); %#ok<*ASGLU>
    epos=pos+ntst;
    %% find other local minima, sorted according to value, if nconnections>1
    % ibounds are the indices at which the orbit is supposed to be cut
    imin=find(diff(sign(diff(npd(pos:epos))))>0)+pos;
    nmin=length(imin);
    assert(nmin>=options.connections-1,'dde_hcli_from_psol:connections',...
        ['dde_hcli_from_psol: requested %d connecting orbits but found only',...
        '%d local minima of derivative'],options.connections,nmin+1);
    [dum,is]=sort(npd(imin));
    imin=imin(is);
    imin=imin(1:options.connections-1);
    imin=sort(imin);
    ibounds=[pos,imin; imin,epos];
else
    %% user has provided hints
    ibounds=[options.eqtime_ind;...
        options.eqtime_ind(2:end),options.eqtime_ind(1)+ntst];
end
for i=size(ibounds,2):-1:1
    %% cut orbit at requested points
    coll=dde_psol_cut(psol,ibounds(:,i));
    hcli(i)=dde_hcli_create('point',coll,...
        'x1',coll.profile(:,1),'x2',coll.profile(:,end));
    %% add linear equilibrium information
    [hcli(i),stst(:,i)]=dde_hcli_from_hcli(hcli(i),pass_on{:});
end
end
%%
function coll=dde_psol_cut(psol,bd)
[psol,submesh]=dde_coll_check(psol);
icoarse=1:psol.degree:length(psol.mesh);
[icoarse,~,ibd]=unique([bd(:)',icoarse]);
tcoarse=psol.mesh(icoarse(ibd(1):ibd(2)));
t=dde_coll_meshfill(tcoarse,1,'grid',submesh);
prof=dde_coll_eva(psol.profile,psol.mesh,t,psol.degree);
tscal=t-t(1);
T=tscal(end);
tscal=tscal/T;
coll=dde_coll_create('mesh',tscal,'parameter',psol.parameter,...
    'period',T*psol.period,'profile',prof,'degree',psol.degree);
end
