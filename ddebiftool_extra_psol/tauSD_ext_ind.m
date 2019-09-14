function xtau_ind=tauSD_ext_ind(ntau)
%% generate array of indices counting additional delays needed for extended SD-DDEs
%
% suppose the original system has the d delays
% tau(1)=0
% tau(j+1)=sys_tau(j,xx(:,1:j),p) (j=1..d-1)
% The extended system has the following combinations:
% arguments are xs(:,1:d,1:d)
% delays are tau(1:d,1:d)
% tau(1,1)=0
% tau(1,j+1)=tau(j+1,1)=sys_tau(j,xs(:,1:j,1),p) (j=1..d-1)
% tau(k+1,j+1)=sys_tau(k,xs(:,1:k,1),p)+sys_tau(j,xs(:,1:k,j+1),p)
% since tau(1,1)=0 and tau(1,j)=tau(j,1) we need only d*(d-1) combinations
%% numbering is row-wise:
% xs(:,1:d,1)=xx(:,1:d)
% xs(:,1,1:d)=xx(:,1:d)
% xs(:,2:d,j)=xx(:,d+(j-1)*(d-1)+(1:d-1))
%% input
% ntau: number of delays in original system
%
%% output
% xtau_ind: xx(:,xtau_ind(i,:)) accesses columns of xx corresponding to
%  [x(t-tau(i,1)),x(t-tau(i,1)-tau(i,1)),x(t-tau(i,1)-tau(i,3)),...]
%
% $Id: tauSD_ext_ind.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
ntau_ext=ntau*(ntau+1);
%% initialise index array into x
xtau_ind=[(1:ntau+1)',reshape(2:ntau_ext+1,ntau,ntau+1)'];
end
