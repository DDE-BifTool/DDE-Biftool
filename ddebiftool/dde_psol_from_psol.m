%% Create initial period doubled orbit near orbit with Flqouet multiplier -1
%%
function [per2,tangent]=dde_psol_from_psol(point,varargin)
%% Inputs
% 
% * |point|: |'psol'| periodic orbits from which one wants to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |funcs| (mandatory): structure with functions provided by user
% * |'radius'|: initial deviation along period-doubling eigenvector
%
% All other name-value pairs are passed on to output branch.
%% Outputs
%
% * |per2|: periodic orbit with initial small period 2 deviation
% * |tangent|: approximate tangent
%
% $Id: dde_psol_from_psol.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
default={'radius',0.01,'method',getfield(df_mthod('psol'),'stability'),'funcs',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% create branch per2 of periodic solutions starting from an
% approximate period doubling num on a branch br of periodic orbits
% obtain eigenvector closest to -1
tangent=repmat(point,0,1);
%% find eigenvector for eigenvalue closest to -1
if ~isfield(point,'stability') || ~isfield(point.stability,'eigenfuncs') ||...
        isempty(point.stability.eigenfuncs)
    if isempty(options.funcs)
        warning('dde_psol_from_psol:arguments',...
            'dde_psol_from_psol: eigenvectors and equations not provided');
        per2=point;
        return
    end
    stability=p_stabil(options.funcs,point,options.method,'closest',-1,...
        'geteigenfuncs',true,'max_number_of_eigenvalues',1,pass_on{:});
else
    stability=point.stability;
end
[~,ipd]=min(abs(stability.mu+1));
if real(stability.mu(ipd))>=0
    warning('dde_psol_from_psol:eigenvalue',...
        'dde_psol_from_psol: eigenvalue closest to -1 %g+1i%g has non-negative real part',...
        real(stability.mu(ipd)),imag(stability.mu(ipd)));
    per2=point;
    return
end
%% use real part of eigenvector to construct small deviation from point
% call it per2. tangent equals eigenvector
deviation=real(stability.eigenfuncs(ipd).profile);
deviation2=[deviation,-deviation(:,2:end)];
deviation2=deviation2/max(abs(deviation2(:)));
per2=dde_psol_create('point',point,...
    'mesh',[point.mesh/2,0.5+point.mesh(2:end)/2],...
    'period',2*point.period,...
    'profile',[point.profile,point.profile(:,2:end)]+options.radius*deviation2);
tangent=dde_psol_create('point',per2,...
    'parameter',0*per2.parameter,...
    'period',0,...
    'profile',deviation2);
end
