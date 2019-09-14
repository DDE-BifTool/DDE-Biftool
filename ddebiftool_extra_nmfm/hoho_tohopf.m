function [hopf,ind]=hoho_tohopf(funcs,point,freqs) 
%% convert Hopf-Hopf point to Hopf bifurcation point
% function hopf_point=hoho_tohopf(funcs,point {,freqs})
% INPUT:
%   funcs problem functions
%	point: double hopf point
%   optional freqs: frequency to be excluded/included
% OUTPUT:
%	hopf_point: Hopf point with single frequency
%
% $Id: hoho_tohopf.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
hopf = point;
hopf.omega = hopf.nvec.omega(1);
ind=1;
hopf.kind = 'hopf';
if isempty(freqs)
    hopf.omega = point.nvec.omega(1);
elseif ischar(freqs)
    ind=str2double(freqs);
    hopf.omega=point.nvec.omega(ind);
else % exclude whichever omega is closer to freqs
    d1=min(abs(point.nvec.omega(1)-freqs));
    d2=min(abs(point.nvec.omega(2)-freqs));
    if d1<d2
        ind=2;
        hopf.omega=point.nvec.omega(2);
    else
        hopf.omega=point.nvec.omega(1);
    end
end
D=ch_matrix(funcs,point.x,point.parameter,1i*hopf.omega);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
hopf.v=E1(:,i2);
if isfield(point,'stability');
    hopf.stability=point.stability;
end
end
