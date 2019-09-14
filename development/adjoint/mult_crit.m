function [eigval,eigprofile]=mult_crit(funcs,point,method,nremove,closest)
%% find Floquet multiplier closest to unit circle and its mode
% function eig_v=mult_crit(period,point)
% INPUT: 
%   funcs problem functions
%   point periodic solution structure
%   method method parameters for computing stability
%   nremove (optional, default=1), number of trivial Floquet multipliers to
%         remove
% OUTPUT:
%   eigval: critical non-trivial eigenvalue (critical=closest to unit
%   circle)
%	eigprofile: profile of eigenvector (mesh equals point.mesh)
%
% $Id: mult_crit.m 19 2014-04-11 14:15:36Z jan.sieber $
%

%% initialise arguments and functions
if nargin<5
    closest=[]; 
end
if nargin<4
    nremove=1;
end
rhomin=method.minimal_modulus; % only find eigenvalues of absolute value>=rhomin
max_number=method.max_number_of_eigenvalues;
col=method.collocation_parameters;
if isfield(method,'delay_accuracy');
    d_ac=method.delay_accuracy;
else
    d_ac=1e-8;
end
[s,S1]=mult_app(funcs,point,rhomin,max_number,col,d_ac);
[dim,ntst]=size(point.profile);
Floquetmodes=reshape(S1(end-dim*ntst+1:end,:),[dim,ntst,length(s)]);
%% remove trivial Floquet multipliers (multiplier closest to unity)
[dum,ix]=sort(abs(s-1)); %#ok<ASGLU>
s(ix(1:nremove))=0;
%% remove eigenvalues with negative imaginary part
imneg=imag(s)<0;
s(imneg)=0;
%% find critical eigenvalue
if isempty(closest)
    [dum,ix]=sort(abs(log(abs(s)))); %#ok<ASGLU>
    eigval=s(ix(1));
    eigprofile=Floquetmodes(:,:,ix(1));
else
    [dum,ix]=sort(log(abs(s-closest))); %#ok<ASGLU>
    eigval=s(ix(1));
    eigprofile=Floquetmodes(:,:,ix(1));
end
