%% residual & jacobian for psol
%
% $Id: dde_psol_jac_res.m 369 2019-08-27 00:07:02Z jansieber $
%%
function [J,res]=dde_psol_jac_res(funcs,pt,free_par,method,varargin)
default={'collocation_parameters',method.collocation_parameters,...
    'phase_condition',method.phase_condition,'rotationcheck',true,...
    'bcmat',eye(size(pt.profile,1)),'pref',[],'phase_tolerance',1e-8};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
pt=dde_coll_check(pt);
[Js,res]=dde_coll_jac(funcs,pt,free_par,...
    'c',options.collocation_parameters,'matrix',method.matrix,...
    'Dtmat',funcs.lhs_matrix(size(pt.profile,1)),pass_on{:});
neqs=length(res);
%% periodicity condition:
[ind,len]=dde_ind_from_point(pt,free_par);
resbc=pt.profile(:,1)-pt.profile(:,end);
if options.rotationcheck
    resbc=mod(resbc+pi,2*pi)-pi;
end
res=[res;options.bcmat*resbc(:)];
n=length(resbc);
Jbc(:,ind.profile(:,1))=eye(n);
Jbc(:,ind.profile(:,end))=-eye(n);
Js.profile=cat(1,Js.profile,options.bcmat*Jbc);
%% phase condition:
if ischar(options.phase_condition) ||options.phase_condition
    if isempty(options.pref)
        pref=pt;
        resph_to_0=true;
    else
        pref=dde_coll_check(options.pref);
        resph_to_0=false;
        if max(abs(diff(pref.profile,[],2)))<options.phase_tolerance
            resph_to_0=true;
        end
    end
    [resph,Wph]=dde_coll_profile_dot(pref,pt,'derivatives',[1,0]);
    resph=resph+0.5*diff(sum(pref.profile(:,[end,1]).^2,1),[],2);
    Jph_dx=sparse(1,size(Js,2));
    Jph_dx(1:ind.profile(end))=pref.profile(:)'*Wph;
    if resph_to_0
        resph=0;
    end
    res=[res;resph];
    Js.profile=[Js.profile;Jph_dx];
end
%% assemble Jacobian
J(:,ind.profile)=sparse(Js.profile);
J(1:neqs,ind.period)=Js.period;
J(1:neqs,ind.parameter)=Js.parameter;
if strcmp(method.matrix,'full')
    J=full(J);
end
end
