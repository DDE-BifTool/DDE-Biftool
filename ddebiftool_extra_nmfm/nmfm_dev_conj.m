function fout=nmfm_dev_conj(fin)
%% create complex conjugate of history function
%
% $Id: nmfm_dev_conj.m 309 2018-10-28 19:02:42Z jansieber $
%% 
fout=nmfm_dev_fun(conj(fin.v),'lambda',conj(fin.lambda),'t',fin.t);
end
