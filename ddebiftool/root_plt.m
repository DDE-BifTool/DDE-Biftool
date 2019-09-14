%% plot roots of linear eigenvalue problem (equilibria)
%%
function root_plt(l0,l1,n1,varargin)
%% function root_plt(l0,l1,n1)
% INPUT:
%	l0 approximated roots 
%	l1 corrected roots
%	n1 Newton convergence of l1 roots
%   'plotaxis' (named optional, default: gca) axis on which to plot
%
% Different colors for unstable (red), stable (green) and critical(black)
%
% $Id: root_plt.m 296 2018-09-24 21:36:56Z jansieber $
%
default={'plotaxis',gca};
options=dde_set_options(default,varargin,'pass_on');
%% uncorrected roots
l0r=l0(real(l0)>0);
l0g=l0(real(l0)<0);
l0k=l0(real(l0)==0);
if isempty(n1),
  n1=ones(size(l1));
end
l1r=l1(real(l1)>0  & n1~=-1);
l1g=l1(real(l1)<0  & n1~=-1);
l1k=l1(real(l1)==0 & n1~=-1);
plot(options.plotaxis,...
    real(l0r),imag(l0r),'rx',...
    real(l0g),imag(l0g),'gx',...
    real(l0k),imag(l0k),'kx',...
    real(l1r),imag(l1r),'r*',...
    real(l1g),imag(l1g),'g*',...
    real(l1k),imag(l1k),'k*');
%% Draw coordinate axes
xlim=get(options.plotaxis,'xlim');
ylim=get(options.plotaxis,'ylim');
dohold=ishold(options.plotaxis); % determine hold status to restore it later
hold(options.plotaxis,'on');
if xlim(1)<0 && xlim(2)>0,
  plot(options.plotaxis,[0 0],ylim,'b-.'); 
end
plot(options.plotaxis,[xlim(1) xlim(2)],[0 0],'b-.');
if ~dohold
    hold(options.plotaxis,'off');
end    
end
