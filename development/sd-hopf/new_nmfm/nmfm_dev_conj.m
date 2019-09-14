function fout=nmfm_dev_conj(fin)
%% create complex conjugate of history function
%
% $Id$
%% 
fout=nmfm_dev_fun(conj(fin.v),'lambda',conj(fin.lambda),'t',fin.t);
end
