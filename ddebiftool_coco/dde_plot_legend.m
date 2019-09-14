function [obj,names]=dde_plot_legend(ax,varargin)
%% Extract display names from objects in plot
% this function either plots the legend , or returns the list of objects
% with unique DisplayName's.
ch=get(ax,'children');
names=arrayfun(@(c)get(c,'DisplayName'),ch,'uniformoutput',false);
[names,ind]=unique(names);
obj=ch(ind);
if length(varargin)<=1 && ...
        (isempty(varargin) || (islogical(varargin{1}) && varargin{1}))
    legend(ax,obj);
else
    legend(ax,obj,varargin{:});
end
end
