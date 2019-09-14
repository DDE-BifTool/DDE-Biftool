function [hoho, success] = p_tohoho(point)
%% Convert to double Hopf point
% function zeho = p_togenh(point)
% INPUT:
%	point: hopf point
% OUTPUT:
%	hoho: uncorrected starting guess for double-Hopf point
%   success: whether conversion was successful
%
% $Id: p_tohoho.m 73 2014-12-31 10:47:51Z jan.sieber $
%
%%

% Set success
success = 1;
hoho = point;

if strcmp(point.kind, 'hopf')
   if ~isfield(point,'stability') || ~isfield(point.stability,'l1') || isempty(point.stability.l1)
      error('P_TOHOHO: point does not contain stability information.');
   end
   hoho.kind = 'hoho';
   hoho.flag = '';
   if ~isfield(hoho,'nmfm')
      hoho.nmfm = [];
   end
   hoho.omega1 = point.omega;
   hoho = rmfield(hoho,'omega');
   [dum,dum,selroot]=nmfm_smrp([],point,[],...
       'remove_omega',true,'threshold',@(imx)imx>1e-8); %#ok<ASGLU>
   if isempty(selroot)
      fprintf('P_TOHOHO: no second imaginary pair!');
      success = 0;
   end
   omega2 = abs(imag(selroot));
   hoho.omega2 = omega2;
else
   fprintf('P_TOHOHO: only hopf points can be converted into double hopf.\n');
   success = 0;
end
end
