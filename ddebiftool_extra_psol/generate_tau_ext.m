function [ext_tau,xtau_ind,relations]=...
    generate_tau_ext(parameter,orig_tau_ind,varargin)
%% generate array of additional delays needed for extended DDE in periodic fold continuation
%
% If the original delays are tau=[0,tau1,tau2,...] then all delays of extended
% system are the combinations
% Tau=repmat(tau,[length(tau),1])+repmat(tau',[1,length(tau)]) 
%
% and the extended delays ext_tau are the upper triangular matrix of 
% Tau(2:end,2:end) (including the main diagonal), stored as
%
% ext_tau=[Tau(2,2:end),Tau(3,3:end),...]
%
% inputs
% parameter: array of original system parameters
% orig_tau_ind: indices of original delays in the system
% ext_tau_base (opt): intended starting index for additional delays in the
%  extended parameter array
%
% outputs
% ext_tau: current values of additional delays
% xtau_ind: xx(:,xtau_ind(i,:)) accesses columns of xx corresponding to
%  [x(t-tau(i)),x(t-tau(i)-tau(1)),x(t-tau(i)-tau(2)),...]
% relations: matrix of size length(ext_tau) x ext_tau_base+length(ext_tau)
%  assuming the additional delays are stored in the parameter
%  array at the positions ext_tau_base+(1:length(ext_tau)), they are defined
%  via sys_cond through relation*parameter=0
%
% $Id: generate_tau_ext.m 309 2018-10-28 19:02:42Z jansieber $
%

%% Process options
default={'ext_tau_base',length(parameter)+2};
options=dde_set_options(default,varargin);
orig_tau=parameter(orig_tau_ind); % value of original delays
ntau=length(orig_tau_ind);        % number of original delays
ntau_ext=ntau*(ntau+1)/2;         % number of additional delays
%% initialise index array into x
xtau_ind=repmat(1:ntau+1,ntau+1,1);
j=0;
for i=1:ntau+1
    for k=i:ntau+1
        j=j+1;
        xtau_ind(i,k)=j;
        xtau_ind(k,i)=j;
    end
end
%% initialise tau_ext values and set up relations
ext_tau=zeros(1,ntau_ext);        % pre-alloc. array of add. delays
relations=zeros(ntau_ext,options.ext_tau_base+ntau_ext);
j=0;
for i=1:ntau
    for k=i:ntau
        j=j+1;
        ext_tau(j)=orig_tau(i)+orig_tau(k);
        relations(j,orig_tau_ind(i))=relations(j,orig_tau_ind(i))+1;
        relations(j,orig_tau_ind(k))=relations(j,orig_tau_ind(k))+1;
        relations(j,j+options.ext_tau_base)=-1;
    end
end
end
