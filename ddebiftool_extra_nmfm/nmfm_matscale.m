function [imatscal,gc]=nmfm_matscale(imat)
%% scale each column of imat by its gcd
% used to reduce computational effort when computing higher-order
% derivatives via polarization identity
%
% for example
% >> imat=[2,3,4;1,3,2]
% imat =
%      2     3     4
%      1     3     2
% >> [imatscal,gc]=nmfm_matscale(imat)
% imatscal =
%      2     1     2
%      1     1     1
% gc =
%      1     3     2
%
% $Id: nmfm_matscale.m 309 2018-10-28 19:02:42Z jansieber $
%%
for i=size(imat,2):-1:1
    [imatscal(:,i),gc(i)]=nmfm_gcd_loc(imat(:,i));
end
end
function [ivscal,gc]=nmfm_gcd_loc(ivec)
gc=ivec(1);
for i=2:length(ivec)
    gc=gcd(gc,ivec(i));
end
ivscal=ivec/gc;
end
