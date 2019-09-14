function [A,B]=dde_psol_monodromy_mixed(J,iper)
%% Floquet problem matrix for mixed (or sparse) problem
%
% * J: output from coll_dde_jac
% * iper(1:2): indices corresponding to x_1(0) and x_n(1) in BVP interval
%
% * A,B: Flopquet multiplier will be given as solution of A v = B mu v
%
%
% $Id: dde_psol_monodromy_mixed.m 305 2018-10-05 20:31:58Z jansieber $
%%
[s1,s2]=size(J);
A=[J;sparse(s2-s1,s2)];
B=sparse(s2,s2);
nfut=s2-iper(2);
A(s1+(1:nfut),iper(2)+(1:nfut))=speye(nfut);
B(s1+(1:nfut),iper(1)+(0:nfut-1))=speye(nfut);
ipast=s1+nfut;
npast=s2-s1-nfut;
A(ipast+(1:npast),iper(2)+(-npast+1:0))=speye(npast);
B(ipast+(1:npast),1:npast)=speye(npast);
end

