function [genh, success] = p_togenh(point)
%% Convert to uncorrected gneeralized Hopf point
% function genh = p_togenh(point)
% INPUT:
%	point: hopf point
% OUTPUT:
%	genh: uncorrected starting guess for generalized hopf point
%   success: whether conversion was successful
%
% $Id: p_togenh.m 66 2014-12-23 23:20:23Z jan.sieber $
%
%% 

% Set success
success = 1;

genh = point;

if strcmp(point.kind, 'hopf')
   genh.kind = 'genh';
   genh.flag = '';
   if ~isfield(genh,'nmfm')
      genh.nmfm = [];
   end
   if ~isfield(genh.nmfm,'L2')
      genh.nmfm.L2 = NaN;
   end
else
   fprintf('P_TOGENH: only hopf points can be converted into generalized hopf.\n');
   success = 0;
end
end
