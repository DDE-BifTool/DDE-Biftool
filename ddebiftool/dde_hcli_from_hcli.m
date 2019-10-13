function [hcli,stst]=dde_hcli_from_hcli(hcli,varargin)
%% recorrect location and linear stability info at equilibria for hcli
%
% $Id$
%%
default={'funcs',[],'stst_correct',true,'method',[],'free_par',[]};
options=dde_set_options(default,varargin,'pass_on');
assert(~isempty(options.funcs),'dde_hcli_from_psol:funcs',...
    'dde_hcli_from_psol: argument ''funcs'' is needed');
if isempty(options.method)
    options.method=df_mthod('stst','cheb');
end
make_stst=@(x)dde_stst_create('x',x,'parameter',hcli.parameter);
correct=@(p)p_correc(options.funcs,p,options.free_par,[],options.method.point,0,p);
stst(1)=make_stst(hcli.x1);
stst(2)=make_stst(hcli.x2);
if options.stst_correct
    [stst(1),suc(1)]=correct(stst(1));
    [stst(2),suc(2)]=correct(stst(2));
    assert(all(suc),'dde_hcli_from_hcli:correction',...
        'dde_hcli_from_hcli: correction failed for equilibrium x%d',...
        find(~suc,1,'first'))
    hcli.x1=stst(1).x;
    hcli.x2=stst(2).x;
end
%% find unstable eigenvalues and right eigenvectors at left end
stst(1).stability=p_stabil(options.funcs,stst(1),options.method.stability);
assert(...
    ~isempty(stst(1).stability.l1) &&...
    max(real(stst(1).stability.l1))>=0,...
    'dde_hcli_from_hcli:stability',...
    'dde_hcli_from_hcli: x1 has no unstable directions');
lambda_v=stst(1).stability.l1(:);
v_sel=real(lambda_v)>0 & imag(lambda_v)>=0;
lambda_v=lambda_v(v_sel);
hcli.lambda_v=lambda_v;
hcli.v=stst(1).stability.v(:,v_sel);
%% find unstable eigenvalues and conjugate left eigenvectors at right end
stst(2).stability=p_stabil(options.funcs,stst(2),options.method.stability);
lambda_w=stst(2).stability.l1(:);
w_sel=real(lambda_w)>0 & imag(lambda_w)>=0;
lambda_w=lambda_w(w_sel);
hcli.lambda_w=lambda_w;
hcli.w=conj(stst(2).stability.w(:,w_sel));
%% estimate coefficients for initial point in unstable subspace
hcli.alpha=hcli.v\(hcli.profile(:,1)-hcli.x1);
hcli.epsilon=sqrt(real(conj(hcli.alpha).'*hcli.alpha));
hcli.alpha=hcli.alpha/hcli.epsilon;
end
