%% test symbolic differentiation for poscontrol
clear
addpath([pwd(),'/../../ddebiftool']);
addpath([pwd(),'/../../ddebiftool_extra_symbolic']);
if sco_isoctave()
    pkg load symbolic
end
%% Set number of delays and create parameter names as strings
syms x x_tau v v_tau % symbols for states and delayed states
parnames={'phi','psi','r','epsilon'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
syms(parnames{:});     % create symbols for parameters
par=sym(parnames);
%% rescaled parameters
alpha=sin(pi/2*phi);
beta=cos(pi/2*phi)*cos(2*pi*psi);
gamma=cos(pi/2*phi)*sin(2*pi*psi);
tau_star=sqrt(8-6*epsilon)/2;
a=1+r^6*alpha;
b=tau_star+r^2*beta*tau_star/3;
tau=b+gamma*tau_star*r^4/3;
%% DDE
dxdt=v;
dvdt=sin(x)-cos(x)*(a*x_tau+b*v_tau);
%% generate code
[fstr,fderivs]=dde_sym2funcs(...
    [dxdt;dvdt],...
    [x,x_tau;v,v_tau],...
    par,...
    'sd_delay',tau,...
    'directional_derivative',true,'filename','sym_balancing');
