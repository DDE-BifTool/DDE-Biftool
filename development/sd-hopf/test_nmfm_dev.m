%% test nmfm call and deriv
clear
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_psol/',...
    '../../ddebiftool_extra_nmfm/',...    
    '../../ddebiftool_utilities/');
addpath('new_nmfm/');
v=reshape(1:6,[2,3]);
lambda=[0,1i,1i];
it=[0,0,1];
f=nmfm_dev_fun(v,'lambda',lambda,'t',[0,0,1]);
theta=linspace(-5,0,200);
x=nmfm_dev_call(f,theta,'c','imag');
dx=nmfm_dev_call(f,theta,'deriv',1,'c','imag');
dt=0.5*(theta(2:end)+theta(1:end-1));
dtx=diff(x,[],2)./[diff(theta,[],2);diff(theta,[],2)];
plot(dt,dtx,theta,dx,'.');