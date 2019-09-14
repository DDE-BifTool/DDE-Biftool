function tau_eva=p_tau(funcs,point,d_nr,t)
%% evaluate (state-dependent) delays along orbit
% function [tau_eva]=p_tau(point,d_nr,t)
% INPUT:
%       point a point
%       d_nr number(s) of delay(s) to evaluate
%       t (optional) point(s) where to evaluate
% OUTPUT:
%       tau_eva value(s) of evaluated delay(s)
%
% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: p_tau.m 19 2014-04-11 14:15:36Z jan.sieber $
%
sys_tau=funcs.sys_tau;

d_nr_eva=length(d_nr); % number of delays to evaluate
max_d_nr=max(d_nr);    % maximum number of delay to evaluate
tau_eva=[];

if strcmp(point.kind,'stst') || strcmp(point.kind,'hopf') || strcmp(point.kind,'fold')
    x=point.x;
    xx=x(:,ones(max_d_nr+1,1));
    tau_eva=NaN(1,d_nr_eva);
    for jj=1:max_d_nr
        tau=sys_tau(jj,xx(:,1:jj),point.parameter);
        for ii=1:d_nr_eva
            if jj==d_nr(ii)
                tau_eva(ii)=tau;
            end
        end
    end
elseif strcmp(point.kind,'psol')
    n=size(point.profile,1);
    mm=point.degree;
    ll=(size(point.profile,2)-1)/mm;
    if isempty(point.mesh)
        mesh=0:1/(ll*mm):1;
    else
        mesh=point.mesh;
    end;
    if nargin>=4
        mesh_eva=t;
        l_me=length(t);
    else
        mesh_eva=mesh;
        l_me=length(mesh);
    end;
    % compute delays
    tau_eva=zeros(1+max_d_nr,l_me);
    xx=NaN(n,l_me,max_d_nr+1);
    % compute delays at mesh_eva
    for jj=1:max_d_nr
        % compute xx at mesh_eva-tau(jj)
        xx(:,:,jj)=psol_eva(point,mesh_eva-tau_eva(jj,:)/point.period);
        for k=1:l_me
            tau_eva(jj+1,k)=sys_tau(jj,reshape(xx(:,k,1:jj),n,jj),point.parameter);
        end
    end
    tau_eva=tau_eva(2:end,:);
    [dum,sel]=ismember(d_nr,1:max_d_nr); %#ok<ASGLU>
    tau_eva=tau_eva(sel,:);
end
end
