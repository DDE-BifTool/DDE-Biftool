function [psolbr,augmented] = nmfm_POfold_from_genh_init(funcs,genh,radius,freepars,varargin)
%% Initialize branch for continuing the Limit point of cycles curve
% emanating from the the generalized-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: nmfm_POfold_from_genh_init.m 373 2019-09-03 22:13:51Z jansieber $
%%
%% create template for periodic solution
psol=p_topsol([],setfield(genh,'kind','hopf'),'radius',radius(1),...
    varargin{:}); %#ok<SFLD>
%% normal form coefficients
n=length(genh.x);
get_n=@(x)x(1:n);
dev_eval=@(x)get_n(nmfm_dev_call(x,0));
phi=dev_eval(genh.nmfm.phi);
h1001=dev_eval(genh.nmfm.h1001);
K01=genh.nmfm.K01;
h0001=dev_eval(genh.nmfm.h0001);
h11=dev_eval(genh.nmfm.h11);
h20=dev_eval(genh.nmfm.h20);
h30=dev_eval(genh.nmfm.h30);
h21=dev_eval(genh.nmfm.h21);
L2=genh.nmfm.L2;
c1=genh.nmfm.c1;
c2=genh.nmfm.c2;
b12=genh.nmfm.b12;
%% parameter expansion
dpar=[genh.parameter(freepars)'];%, zeros(length(freepars),1), -2*real(c2)*K01];
%% frequency expansion
d_om=[genh.omega];%, 0, imag(c1)-2*real(c2)*b12];
% %% profile expansion
% zn=zeros(n,1);
% profile={...
%     [genh.x, zn,    h11-2*L2*h0001],...
%     [zn,     2*phi, zn,            -4*L2*h1001+h21],...
%     [zn,     zn,    h20],...
%     [zn,     zn,    zn,             1/3*h30]};
%% profile expansion
zn=zeros(n,1);
profile={...
    [genh.x, zn,    h11-2*real(c2)*h0001],...
    [zn,     2*phi, zn],...
    [zn,     zn,    h20]};
%% approximation to limit cycle
psolbr=nmfm_psol_from_C2(funcs,profile,dpar,d_om,freepars,radius,psol);
augmented=false;
end
