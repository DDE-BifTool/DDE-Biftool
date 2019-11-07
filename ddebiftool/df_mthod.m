function method=df_mthod(varargin)
%% create methed structure with default parameters for continuation, solving and stability
% function method=df_mthod(kind,flag_newhheur)
% INPUT:
%	kind the kind of default method wanted
%   opt: discretization (optional, default: 2) boolean: use Chebyshev approx
%                     Only used if kind==stst, alternatives: 0=bdf, 1=mxo
%                     (Verheyden)
% OUTPUT:
%	method default method
% COMMENT:
%       The order of the LMS method used in the computation of the
%       characteristic roots is hard coded in this file.
%  
% Update on 05/03/2007 ("flag_newhheur" <=> (imag(method.stability.lms_parameter_rho)~=0) )
%
% $Id: df_mthod.m 369 2019-08-27 00:07:02Z jansieber $
%
%% correction of argument list for compatibility
% (formerly first arg funcs no longer needed)
flag_newhheur=2; % compatibility to Verheyden format (but switching away from Verheyden's method)
if length(varargin)==1
    kind=varargin{1};
elseif length(varargin)==2 && isstruct(varargin{1})
    kind=varargin{2};
elseif length(varargin)==2 && ischar(varargin{1})
    kind=varargin{1};
    flag_newhheur=varargin{2};
elseif length(varargin)>=3
    kind=varargin{2};
    flag_newhheur=varargin{3};
end
try
    method=dde_apply({'dde_',kind,'_method'},varargin{3:end});
    return
catch ME
    if ~strcmp(ME.identifier,'dde_apply:function')
        rethrow(ME);
    end
end
%% continuation parameters
method.continuation.steplength_condition=1; % use steplength condition
method.continuation.plot=1; % plot new points
method.continuation.prediction=1; % use secant prediction
method.continuation.steplength_growth_factor=1.2; % grow steplength with factor 1.2
method.continuation.steplength_minimum=0; % minimal stepsize
method.continuation.plot_progress=1; % plot progress gradually
method.continuation.plot_measure=[]; % use default plot measures
method.continuation.halt_before_reject=0; % rejection of points allowed
method.continuation.minimal_angle=-1; % minimal cos(angle) between prediction and correction

%% Bifurcation detection
method.bifurcation.radial_tolerance_factor = 0.25;
method.bifurcation.minimal_real_part = -0.1;
method.bifurcation.correction_tolerance = 1e-6;
method.bifurcation.secant_iterations = 30;
method.bifurcation.secant_tolerance = 1e-6;
method.bifurcation.print = 0;
method.bifurcation.imagthreshold = 1e-6;
method.bifurcation.monitor_eigenvalues = 0;
method.bifurcation.plot_testfunctions = 0;
method.bifurcation.pause_on_bifurcation= 0;
%% Newton iteration
method.point.newton_max_iterations=5;
if strcmp(kind,'hcli')
    method.point.newton_max_iterations=10;
end
method.point.newton_nmon_iterations=1;
method.point.preprocess='';
method.point.postprocess='';
method.point.extra_condition=0;
method.point.print_residual_info=0;
method.point.jacobian_nonsquare=false;
method.point.remesh=false;
method.point.delay_zero_prep='stst';
switch kind
    case 'stst'
        % stst
        method.point.halting_accuracy=1e-10;
        method.point.minimal_accuracy=1e-8;
    case 'fold'
        method.point.halting_accuracy=1e-9;
        method.point.minimal_accuracy=1e-7;
    case 'hopf'
        method.point.preprocess='dde_jac2square_preprocess';
        method.point.postprocess='dde_hopf_postprocess';
        method.point.halting_accuracy=1e-9;
        method.point.minimal_accuracy=1e-7;
    case {'psol','hcli'}
        % point
        method.point.halting_accuracy=1e-8;
        method.point.minimal_accuracy=1e-6;
        method.point.phase_condition=1;
        method.point.collocation_parameters=[];
        method.point.adapt_mesh_before_correct=0;
        method.point.adapt_mesh_after_correct=3;
        method.point.remesh=true;
        method.point.delay_zero_prep='coll';
        method.point.matrix='full';
    otherwise
        [dum1,dum2,kindparent]=feval(['dde_',kind,'_create']); %#ok<ASGLU>
        if ~isempty(kindparent)
            method=df_mthod(kindparent,flag_newhheur,varargin{4:end});
            return
        else
            warning('df_mthod: kind %s not recognized',kind);
        end
end
%% Stability/spectrum
discretizations={'bdf','mxo','cheb'};
method.stability.delay_accuracy=-1e-8;
method.stability.max_number_of_eigenvalues=20;
method.stability.root_accuracy=1e-6;
method.stability.fill=[];
switch kind
    case {'stst','hopf','fold'}
        if isnumeric(flag_newhheur)
            discretization=discretizations{flag_newhheur+1};
        elseif ischar(flag_newhheur)
            discretization=flag_newhheur;
        end
        method.stability.discretization=discretization;
        method.stability.delay_accuracy=-1e-8;
        method.stability.root_accuracy=1e-6;
        switch discretization
            case{'bdf','mxo'}
                order=4;
                delta_region=0.1;
                [alpha,beta]=dde_stst_time_lms(discretization,4);
                method.stability.lms_parameter_alpha=alpha;
                method.stability.lms_parameter_beta=beta;
                switch discretization
                    case 'mxo'
                        method.stability.lms_parameter_rho=dde_stst_mxo_ellipse(order,delta_region);
                        method.stability.newheuristics_tests=2000;
                    case 'bdf'
                        method.stability.lms_parameter_rho=dde_stst_bdf_safety(alpha,beta,0.01,0.01);
                end
                method.stability.interpolation_order=order;
                method.stability.minimal_time_step=0.01;
                method.stability.maximal_time_step=0.1;
                method.stability.max_newton_iterations=6;
                method.stability.remove_unconverged_roots=1;
                method.stability.minimal_real_part=[];
            case 'cheb'
                method.stability.max_number_of_eigenvalues=20;
                method.stability.min_number_of_eigenvalues=[];
                method.stability.ncheb=[]; % initial 
                method.stability.inisize=[];
                method.stability.maxsize=200;
                method.stability.minimal_real_part=-1;       
                method.stability.scale_w=true;
                method.stability.nearest=[];
                method.stability.tauscal_limit=1e-5;
                method.stability.discard_accuracy_factor=1e5; % if finite, roots with err>root_accuracy*this will be discarded
        end
    case 'psol'   % multipliers
        method.stability.eigmatrix='full';
        method.stability.closest=[];
        method.stability.fill=[];
        method.stability.minimal_modulus=0;
        method.stability.root_accuracy=1e-6;
        method.stability.collocation_parameters=[];
        method.stability.geteigenfuncs=false;
    case 'hcli'
        method=rmfield(method,'stability');
end

