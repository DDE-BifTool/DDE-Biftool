function cusp = nmfm_cusp(funcs, cusp, varargin)
%% Cusp normal form coefficients
%
% This function assumes that the input cusp is at least of type fold.
%
% $Id: nmfm_cusp.m 309 2018-10-28 19:02:42Z jansieber $
%% 
default={'nullpoint',[],'free_pars',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[p,q]=nmfm_nullvector(funcs,cusp,0,'nullpoint',options.nullpoint,pass_on{:});

Delta = ch_matrix(funcs,cusp.x,cusp.parameter,0);
DDelta = ch_matrix(funcs,cusp.x,cusp.parameter,0,'deri',1);
par0=cusp.parameter(:)*0;
dev0=@(v)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))]);
F=nmfm_deriv_define(funcs,cusp,...
    'free_pars',options.free_pars,pass_on{:});
%% eigenfunctions
phi = dev0(q);

Bord = [Delta q; p 0];
B_phi2=F.B(phi,phi);
DinvD2f = Bord\[B_phi2; 0];
DinvD2f= DinvD2f(1:end-1);

h2 = dev0(-DinvD2f + p*DDelta*DinvD2f*q);

%% critical normal forms coefficients
b=1/2*p*B_phi2;
c=1/6*p*(3*F.B(phi,h2)+F.C(phi,phi,phi));

% b=b/alpha;
% c=c/alpha;

cusp.nmfm.b=b;
cusp.nmfm.c=c;
cusp.nvec.q=q;
cusp.nvec.p=p;


end
