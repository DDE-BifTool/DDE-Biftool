function [psolbr,augmented] = nmfm_TorusBifurcation_from_hoho_init(funcs,hoho,radius,freepars,varargin)
%% Initialize branch for continuing the Neimark-Sacker curves
% emanating from the the Hopf-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% (adapted by JS)
%
% $Id: nmfm_TorusBifurcation_from_hoho_init.m 369 2019-08-27 00:07:02Z jansieber $
%
%%
%% create template for periodic solution
psoltemplate=p_topsol([],setfield(hoho,'kind','hopf'),'radius',radius(1),...
    varargin{:}); %#ok<SFLD>
%% normal form coefficients & eigenvectors
n=length(hoho.x);
get_n=@(x)x(1:n);
dev_eval=@(x)get_n(nmfm_dev_call(x,0));
beta=-[...
    hoho.nmfm.g2100, hoho.nmfm.g1011;...
    hoho.nmfm.g1110, hoho.nmfm.g0021];
rebeta=real(beta);
imbeta=imag(beta);
h0011=dev_eval(hoho.nmfm.h0011);
delta=  [dev_eval(hoho.nmfm.h2000),    dev_eval(hoho.nmfm.h0020)];
h0000mu=[dev_eval(hoho.nmfm.h0000(1)), dev_eval(hoho.nmfm.h0000(2))];
omega=hoho.nvec.omega;
b=hoho.nmfm.b;
K=hoho.nmfm.K;
q=hoho.nvec.q;
zn=zeros(n,1);
rotation_numbers=exp(1i*abs(mod(2*omega(2:-1:1)./omega(1:2)+1,2)-1)*pi);
for ib=2:-1:1
    %% parameter expansion
    dpar=[hoho.parameter(freepars)',zeros(length(freepars),1),K*rebeta(:,ib)];
    %% frequency expansion
    d_om=[omega(ib),0,(b(:,ib)'*rebeta(:,ib)-imbeta(ib,ib))];
    %% profile expansion
    profile={...
        [hoho.x, zn,        h0011+h0000mu*rebeta(:,ib)],...
        [zn,     2*q(:,ib)],...
        [zn,     zn,        delta(:,ib)]};
    psolbr(ib)=nmfm_psol_from_C2(funcs,profile,dpar,d_om,freepars,radius,psoltemplate);
    psolbr(ib)=replace_branch_pars(psolbr(ib),psolbr(ib).parameter.free,varargin);
end
augmented={...
    {'closest',rotation_numbers(1)},...
    {'closest',rotation_numbers(2)}};
end
