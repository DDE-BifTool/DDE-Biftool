function [hoho, success] = p_tohoho(point)
%% Convert to double Hopf point
% function hoho = p_tohoho(point)
% INPUT:
%	point: hopf point
% OUTPUT:
%	hoho: uncorrected starting guess for double-Hopf point
%   success: whether conversion was successful
%
% $Id: p_tohoho.m 314 2019-01-24 14:28:23Z mmbosschaert $
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
      hoho.nmfm = struct();
   end
   if ~isfield(hoho,'nvec') || isempty(hoho.nvec) || ~isfield(hoho.nvec,'omega')
       hoho.nvec.omega=point.omega;
   end
   if length(hoho.nvec.omega)==2
       return
   end
   [dum,dum,selroot]=nmfm_smrp([],point,[],...
       'remove_omega',true,'threshold',@(imx)imx>1e-8); %#ok<ASGLU>
   if isempty(selroot)
      fprintf('P_TOHOHO: no second imaginary pair!');
      success = 0;
   end
   hoho.nvec.omega(2) = abs(imag(selroot));
else
   fprintf('P_TOHOHO: only hopf points can be converted into double hopf.\n');
   success = 0;
end

% switch order of omega's if needed
if hoho.nvec.omega(1) > hoho.nvec.omega(2)
    hoho.nvec.omega(:)=hoho.nvec.omega([2 1]);
end
