function hopf=zeho_tohopf(funcs,point,freqs) %#ok<INUSD>
%% convert zero-Hopf point to Hopf bifurcation point
% function hopf_point=zeho_tohopf(funcs,point {,freqs})
% INPUT:
%   funcs problem functions
%	point with stability information 
%   optional freqs: frequency to be excluded from consideration
% OUTPUT:
%	hopf_point: hopf point
%
% $Id: zeho_tohopf.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
hopf = point;
hopf.kind='hopf';
x = point.x;
D=ch_matrix(funcs,x,point.parameter,1i*point.nvec.omega);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
hopf.v=E1(:,i2);
hopf.omega=point.nvec.omega;
if isfield(point,'stability');
    hopf.stability=point.stability;
end
end
