function [genh, success] = p_togenh(point)
%% Convert to uncorrected gneeralized Hopf point
% function genh = p_togenh(point)
% INPUT:
%	point: hopf point
% OUTPUT:
%	genh: uncorrected starting guess for generalized hopf point
%   success: whether conversion was successful
%
% $Id: p_togenh.m 309 2018-10-28 19:02:42Z jansieber $
%
%% 

% Set success
success = 1;

genh = point;

if strcmp(point.kind, 'hopf')
   genh.kind = 'genh';
   genh.flag = '';
   if ~isfield(genh,'nmfm')
      genh.nmfm = struct();
   end
   if ~isfield(genh.nmfm,'L2')
      genh.nmfm.L2 = NaN;
   end
else
   fprintf('P_TOGENH: only hopf points can be converted into generalized hopf.\n');
   success = 0;
end
end
