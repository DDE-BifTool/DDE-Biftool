%% test composite derivative
%s=load('poscontrol_results.mat','hopf_wbifs','hopf');
s=load('nested_single_results.mat','hopfs');
hopf_wbifs=s.hopfs;
% ntau=3;
% parnames={'tau0','s0','c','K'};
% cind=[parnames;num2cell(1:length(parnames))];
% ind=struct(cind{:});
%rhs=@(x,s,tau0,s0,c,K)[-K*c/2*(s(2,:)-s0);...
%    (-K*c/2*(s(4,:)+s(2,:)-2*s0)-c*s(1,:)+x(3,:)+x(1,:))/c+K*c/2*(s(4,:)-s0)];
%rhs_xp=@(xx,par)rhs(reshape(xx(1,:,:),ntau+1,[]),reshape(xx(2,:,:),ntau+1,[]),par(1),par(2),par(3),par(4));
funcs=set_symfuncs('sym_nested');
nfuncs=set_funcs('sys_rhs',@(xx,p)-xx(:,2,:),'sys_tau',@(it,xx,p)p+xx(:,1,:),'sys_ntau',@()1,...
    'x_vectorized',true);
%%
for i=1:length(s.hopfs.point)
    hopf=s.hopfs.point(i);
    newpoint=nmfm_hopf(funcs,hopf);
    n2point=nmfm_hopf(nfuncs,hopf);
    fprintf('L1 num=%g, L1 nest=%g, diff=%g\n',...
        n2point.nmfm.L1,newpoint.nmfm.L1,n2point.nmfm.L1-newpoint.nmfm.L1)
end