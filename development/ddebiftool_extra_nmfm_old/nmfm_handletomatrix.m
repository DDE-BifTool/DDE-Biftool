function PHI = nmfm_handletomatrix(fn, arg)
%% Call single-argument function with array of arguments
% INPUT:
%   fn: function handle
%   arg: argument vector for function
% OUTPUT:
%   PHI: n by r matrix
%
% $Id: nmfm_handletomatrix.m 65 2014-12-23 22:37:36Z jan.sieber $
%
%%
for k = length(arg):-1:1
    phivec = fn(arg(k));
    PHI(:,k) = phivec(:);
end
end

