function [mats,factors]=nmfm_matcollect(mats,factors)
%% collect equal polarization combinations
%
% [mats,factors]=nmfm_matcollect(mats,factors)
%
% to reduce computational effort if polarization identity creates many
% (equal) terms. For example, if mats is 
% [ 1,1,1;
%   0,0,1] and factors is 
% [ 2,3,4] then the two [1;0] columns of mats can be combined such that 
% output will be mats=
% [ 1,1;
%   0,1] and factors [5,4].
%
% $Id$
%% 
[mats,~,ind]=unique(mats.','rows');
mats=mats';
factors=accumarray(ind,factors).';
end