function theme = dde_plot_theme(branchtype,prefix)
%DDE_PLOT_THEME    Default plot theme, adjusted from for COCO's PO toolbox.
%
%   DDE_PLOT_THEME with no arguments displays a list of default themes for
%   all branch types
%
%   DDE_PLOT_THEME(branchtype) returns a struct defining the plot theme for
%   a branchtype. This theme struct contains only settings that
%   override or augment the default COCO plot theme. The full default theme
%   struct is computed in COCO_PLOT_THEME as
%   COCO_MERGE(COCO_PLOT_THEME,DDE_PLOT_THEME(branchtype)).

% $Id: dde_plot_theme.m 346 2019-05-13 05:41:50Z jansieber $

if nargin<1
  display_themes();
else
  if nargin<2
      prefix='';
  end
  theme = struct();
  theme.sol.RO   = {'k-', 'LineWidth', 1, 'Marker', '.'};
  theme.plot_sol = @po_plot_sol;
  
  switch branchtype
    
      case {'dde.stst','stst'}
          theme.bd.col2  = coco_get_id(prefix,'xnorm');
          lspec_s = {'b-', 'LineWidth', 1,'DisplayName','stst stable'};
          lspec_u = {'b--', 'Linewidth', 1,'DisplayName','stst unstable'};
          theme.lspec    = {lspec_s, lspec_u};
          theme.ustab    = coco_get_id(prefix,'USTAB');
          theme.ustabfun = @(x) (x>=1)+1;
          theme.special  = {'HB', 'SN', 'BP'};
          theme.usept    = {'HB', 'SN', 'FP', 'BP'};
          theme.HB       = {'ko', 'MarkerFaceColor', 'c', 'MarkerSize', 8,'DisplayName','HB stst'};
          theme.SN       = {'kd', 'MarkerFaceColor', 'g', 'MarkerSize', 8,'DisplayName','SN stst'};
          theme.BP       = {'kd', 'MarkerFaceColor', 'w', 'MarkerSize', 8,'DisplayName','BP stst'};
          theme.NSA      = {'kd', 'MarkerFaceColor', 'm', 'MarkerSize', 8,'DisplayName','NSA stst'};
          theme.sol.HB   = {'mo', 'LineWidth', 2,'DisplayName','HB stst'};
          theme.sol.SN   = {'g-', 'LineWidth', 2,'DisplayName','SN stst'};

      case {'dde.psol','psol'}
          theme.bd.col2  = coco_get_id(prefix,'xnorm');
          lspec_s = {'r-', 'LineWidth', 2,'DisplayName','psol stable'};
          lspec_u = {'r--', 'Linewidth', 2,'DisplayName','psol unstable'};
          theme.lspec    = {lspec_s, lspec_u};
          theme.ustab    = coco_get_id(prefix,'USTAB');
          theme.ustabfun = @(x) (x>=1)+1;
          theme.usept    = {'SN', 'PD', 'TR', 'FP', 'BP'};
          theme.special    =  {'SN', 'PD', 'TR', 'BP'};
          theme.SN       = {'kd', 'MarkerFaceColor', 'g', 'MarkerSize', 8,'DisplayName','SN psol'};
          theme.PD       = {'kd', 'MarkerFaceColor', 'r', 'MarkerSize', 8,'DisplayName','PD psol'};
          theme.TR       = {'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8,'DisplayName','TR psol'};
          theme.BP       = {'kd', 'MarkerFaceColor', 'w', 'MarkerSize', 8,'DisplayName','BP psol'};
          theme.NSA      = {'kd', 'MarkerFaceColor', 'm', 'MarkerSize', 8,'DisplayName','NSA psol'};
          theme.sol.SN   = {'g-', 'LineWidth', 2,'DisplayName','SN psol'};
          theme.sol.PD   = {'r-', 'LineWidth', 2,'DisplayName','PD psol'};
          theme.sol.TR   = {'b-', 'LineWidth', 2,'DisplayName','TR psol'};
          
      case {'dde.hopf','hopf'}
          theme.lspec   = {'k-', 'LineWidth', 2,'DisplayName','Hopf'};
          theme.usept   = {'SN', 'HB', 'BP'};
          theme.special = theme.usept;
          theme.HB       = {'ko', 'MarkerFaceColor', 'c', 'MarkerSize', 8,'DisplayName','Hopf-Hopf'};
          theme.SN       = {'kd', 'MarkerFaceColor', 'g', 'MarkerSize', 8,'DisplayName','Hopf-Fold'};
          
      case {'dde.fold','fold'}
          theme.lspec   = {'g-', 'LineWidth', 2,'DisplayName','Fold stst'};
          theme.usept   = {'HB', 'BP'};
          theme.special = theme.usept;
          theme.HB       = {'kd', 'MarkerFaceColor', 'g', 'MarkerSize', 8,'DisplayName','Hopf-Fold'};

      case {'dde.POfold','POfold'}
          theme.lspec   = {'-', 'LineWidth', 2,'color',[0.5,0.5,0],'DisplayName','SN psol'};
          
      case {'dde.PD','PD'}
          theme.lspec   = {'m-', 'LineWidth', 2,'DisplayName','PD'};
          
      case {'dde.torus','torus'}
          theme.lspec   = {'c-', 'LineWidth', 2,'DisplayName','TR'};
          
      otherwise
          error('%s: unknown solution branch type ''%s''', mfilename, branchtype);
  end
  
end
end

function display_theme(bt)
tb_info.tb = 'dde';
tb_info.dde.branch_type = coco_get_id(tb_info.tb,bt);
btname = coco_plot_theme(tb_info);
disp(' ');
disp(['',tb_info.dde.branch_type,''' =']);
disp(' ');
disp(btname);
end
function display_themes()
%tb_info.tb = 'dde';

fprintf(' default plotting themes:\n');

display_theme('stst');
display_theme('psol');
display_theme('fold');
display_theme('hopf');
display_theme('PD');
display_theme('torus');

% 
% tb_info.dde.branch_type = 'dde.HB';
% HB = coco_plot_theme(tb_info);
% disp(' ');
% disp('''dde.HB'' =');
% disp(' ');
% disp(HB);
% 
% tb_info.dde.branch_type = 'dde.PD';
% PD = coco_plot_theme(tb_info);
% disp(' ');
% disp('''dde.PD'' =');
% disp(' ');
% disp(PD);
% 
% tb_info.dde.branch_type = 'dde.TR';
% TR = coco_plot_theme(tb_info);
% disp(' ');
% disp('''dde.TR'' =');
% disp(' ');
% disp(TR);

end
