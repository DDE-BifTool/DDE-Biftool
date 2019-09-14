function data=dde_psol2psol(funcs,info,varargin)
default={'radius',0.01,'method',info.method.stability,'biftype','psol'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
point=info.point(1);
info=replace_branch_pars(info,[],pass_on);
[info.point,info.tangent]=feval(['dde_init_',options.biftype],funcs,...
    point,info,options.radius,pass_on{:});
data.info=info;
data.funcs=funcs;
end
%%
function [pt,tangent]=dde_init_psol(~,point,info,~,varargin)
pt=point;
tangent=info.tangent;
end
%%
function [per2,tangent]=dde_init_PD(funcs,point,info,radius,varargin)
method=info.method;
stability=p_stabil(funcs,point,method.stability,'closest',-1,...
    'geteigenfuncs',true,'max_number_of_eigenvalues',1,varargin{:});
[~,ipd]=min(abs(stability.mu+1));
if real(stability.mu(ipd))>=0
    error('dde_init_PD:eigenvalue',...
        'dde_init_PD: eigenvalue closest to -1 %g+1i%g has non-negative real part',...
        real(stability.mu(ipd)),imag(stability.mu(ipd)));
end
deviation=real(stability.eigenfuncs(ipd).profile);
deviation2=[deviation,-deviation(:,2:end)];
deviation2=deviation2/max(abs(deviation2(:)));
per2=dde_psol_create('point',point,...
    'mesh',[point.mesh/2,0.5+point.mesh(2:end)/2],...
    'period',2*point.period,...
    'profile',[point.profile,point.profile(:,2:end)]+radius*deviation2);
tangent=dde_psol_create('point',per2,...
    'parameter',0*per2.parameter,...
    'period',0,...
    'profile',deviation2);
end
