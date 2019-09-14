%% SetupMWFold - Initialize continuation of folds of modulated waves
%%
function [pfuncs,pbranch,suc]=SetupMWFold(funcs,branch,ind,varargin)
%% Inputs
% 
% * |funcs|: functions used for rotationally symmetric DDE
% * |branch|: branch of psols along which fold was discovered
% * |ind| number of point close to fold
%
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%% Optional inputs
% 
% * |contpar| (integers default |[]|): index of continuation parameters 
%   (replacing free parameters in argument branch)
% * |hbif| (default |1e-3|): used for finite differencing when approximating
%   linearized system, replace by |funcs.sys_deri| if |funcs.sys_deri| is
%   analytical
% * |correc| (logical, default |true|): apply |p_correc| to first points on fold
%   branch
% * |dir| (integer, default |[]|): which parameter to vary initially along fold
%   branch (|pbranch| has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if dir is non-empty
%
% all other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter
%
% $Id: SetupMWFold.m 374 2019-09-14 14:02:58Z jansieber $
%

%% process options
branch.point=branch.point(ind);
[funcs,branch]=PsolFromPsolbif(funcs,branch);
default={'nullparind',branch.parameter.free(end),...
    'extra_cond',{@(p,pref)sys_cond_MWFold(p,pref,funcs)}};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[pfuncs,pbranch,suc]=SetupPOfold(funcs,branch,1,...
    'nullparind',options.nullparind,...
    'extra_cond',options.extra_cond,pass_on{:});
end
% 
% branch.point=branch.point(ind); % remove all points but approx fold
% % initialize branch of folds (pbranch)
% pbranch=replace_branch_pars(branch,options.contpar,pass_on);
% %% set up numbering and values of additional parameters and artificial delays
% point=dde_psol_create('point',pbranch.point);
% ip.dim=size(point.profile,1);    % dimension of original problem
% ip.omega=length(point.parameter); % number of original system parameters
% ip.beta=ip.omega+1;              % location of add. parameter beta
% ip.period=ip.omega+2;            % location of copy of period
% ip.rho=ip.omega+3;
% ip.null_par=[ip.beta,ip.rho];
% ip.orig_tau=funcs.sys_tau(); % system delays
% % additional delays needed for extended system:
% ip.ext_tau=ip.rho+(1:length(ip.orig_tau));
% %% set up functions of extended system
% % (constant delays only for now)
% user_lhs_num=funcs.lhs_matrix(ip.dim);
% lhs_num=kron(eye(2),user_lhs_num);
% pfuncs=funcs;
% pfuncs.lhs_matrix=@(n)lhs_num;
% pfuncs.get_comp=@(p,component)extract_from_MWFold(p,component,ip.omega);
% pfuncs.sys_rhs=@(x,p)sys_rhs_MWFold(x,p,ip,funcs);
% pfuncs.orig_cond=funcs.sys_cond;
% pfuncs.orig_cond_reference=funcs.sys_cond_reference;
% pfuncs.sys_tau=@()[ip.orig_tau,ip.ext_tau];
% pfuncs.delayed_derivs=[zeros(1,length(ip.orig_tau)),ones(1,length(ip.ext_tau))];
% pfuncs=pfuncs.add_deriv(pfuncs,@(x,p,dx,dp)pfuncs.sys_rhs(x,p),0);
% pfuncs.kind='POfold';
% pfuncs.userfuncs=funcs;
% %% required amendments of structures for extended system
% pbranch.parameter.free=[pbranch.parameter.free,...
%     ip.beta,ip.period,ip.rho,ip.ext_tau];
% pbranch.method.point.extra_condition=1;
% %% create initial guess for correction
% pfoldini0=MWFoldInit(funcs,point,branch.method.point,branch.parameter.free);
% %% correct initial guess and find 2nd point along branch if desired
% [pbranch,suc]=correct_ini(pfuncs,pbranch,pfoldini0,...
%     options.dir,options.step,options.correc);
% end
% 
% %% crude initial guess
% function pfoldini=MWFoldInit(funcs,point,method,free_par_ind)
% pfoldini=point;
% J=p_correc_rhs(funcs,method,point,free_par_ind,'pref',point,'output','J');
% nullvecs=dde_nullspaces_lr(J,'nulldim',1);
% nx=numel(point.profile);
% v=nullvecs(1:nx);
% beta=nullvecs(nx+1);
% rho=nullvecs(end);
% vpoint=p_axpy(0,point,[]);
% vpoint.profile=reshape(v,size(point.profile));
% vpoint.parameter=[beta,rho];
% normv=sqrt(p_dot(vpoint,vpoint,'free_par_ind',[1,2]));
% beta=beta/normv;
% rho=rho/normv;
% vpoint.profile=vpoint.profile/normv;
% pfoldini.profile=[pfoldini.profile;vpoint.profile];
% ext_tau=pfoldini.parameter(funcs.sys_tau());
% pfoldini.parameter=[pfoldini.parameter,beta,pfoldini.period,rho,ext_tau];
% end
% %% extract components from MWFold
% function result_array=extract_from_MWFold(pfold_array,component,npar)
% %% check if input is branch rather than point array
% if ~isfield(pfold_array,'kind') && isfield(pfold_array,'point')
%     pfold_array=pfold_array.point;
% end
% dim=size(pfold_array(1).profile,1)/2;
% for i=1:length(pfold_array)
%     pfold=pfold_array(i);
%     switch component
%         case 'kind'
%             result='MWFold';
%         case {'solution','solution_for_stability'}
%             result=pfold;
%             result.profile=result.profile(1:dim,:);
%             result.parameter=result.parameter(1:npar);
%         case 'nullvector'
%             result=pfold;
%             result.profile=result.profile(dim+1:end,:);
%             result.parameter=result.parameter(npar+1:3);
%         case 'delays'
%             result.values=pfold.parameter(npar+4:end);
%         case 'omega'
%             result=pfold.parameter(npar);
%         case 'beta'
%             result=pfold.parameter(npar+1);
%         case 'rho'
%             result=pfold.parameter(npar+3);
%         otherwise
%             error('MWFold:unknown','component %s unknown',component);
%     end
%     result_array(i)=result; %#ok<AGROW>
% end
% end
