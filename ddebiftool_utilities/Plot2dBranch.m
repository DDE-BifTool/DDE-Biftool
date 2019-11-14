%% Plot branch solution branch in bifurcation diagram
%
%% Immportant optional inputs
%
% * |'name'| label used for branch in legend
% * |'ax'| axis to which plot is added (default gca)
% * |'color'| base color for branch (will be gradually paler for unstabile
% eigenvalues) ,default is from line palette
% * |'funcs'|: problem functions, needed if stability not yet computed
% * |'parameter'|: parameters to be used for plotting (default is
% |branch.parameter.free|)
% * |'x'|, |'y'| (overriding |'parameter'|): functions @(p) -> R, to
% determine coordinates from point p
%
%% Outputs
%
% * |lg|: legend handle and text, as filled up to now. To be passed on to
% next call
% * |ax|: axes handle where plot was put
%
% $Id: Plot2dBranch.m 369 2019-08-27 00:07:02Z jansieber $
function [lg,ax]=Plot2dBranch(branch,varargin)
%% 
%#ok<*AGROW>
%% Call recursively for multiple branches
if iscell(branch) && length(branch)>=1
    [lg,ax]=Plot2dBranch(branch{1},varargin{:});
    ish=ishold(ax);
    hold(ax,'on');
    for i=2:length(branch)
        [lg,ax]=Plot2dBranch(branch{i},varargin{:});
    end
    if ~ish
        hold(ax,'off');
    end
    return
elseif isempty(branch)
    lg={};
    ax=[];
    return
end
default={'name',[],'lgname',[],'ax',gca,'funcs',[],'color',[],'linewidth',2,'markersize',8,'bifcolor','k',...
    'stability',0.5,'parameter',branch.parameter.free,'markersymbols','x*sdph+o',...
    'x',[],'y',[],'oldlegend',{[],{}},'pointtype_list',@pointtype_list,...
    'max_nunst',Inf};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if options.stability && isfield(branch.method,'stability')
    [nunst,~,~,pts]=GetStability(branch,'funcs',options.funcs,...
        'exclude_trivial',true,'pointtype_list',options.pointtype_list,pass_on{:});
    stprint=@(name,i)stprint_loc(name,i,options.max_nunst);
else
    pts=branch.point;
    nunst=zeros(size(pts));
    stprint=@(name,i)stprint_loc(name,i,0);
end
nunst=min(nunst,options.max_nunst);
nunst_unique=unique(nunst);
par=@(sel,p_ind)arrayfun(@(x)x.parameter(p_ind),pts(sel));
if isempty(options.x)
    p1=@(pt)par(pt,options.parameter(1));
else
    p1=@(ind)arrayfun(options.x,pts(ind));
end
if isempty(options.y) && length(branch.parameter.free)>=2
    p2=@(pt)par(pt,options.parameter(2));
elseif isempty(options.y) && length(branch.parameter.free)<2
    if isfield(pts(1),'profile')
        p2=@(ind)arrayfun(@(p)max(p.profile(1,:)),pts(ind));
    elseif isfield(pts(1),'x')
        p2=@(ind)arrayfun(@(p)p.x(1),pts(ind));
    else
        error('Plot2dBranch:yaxis','Plot2dBranch: no y axis chosen');
    end
else
    p2=@(ind)arrayfun(options.y,pts(ind));
end    
%% assign colors to curves
typenames=fieldnames(getfield(options.pointtype_list(),'codim')); %#ok<GFLD>
if ~isempty(options.name)
    name=options.name;
elseif ~isempty(options.funcs) && isfield(options.funcs,'kind')
    name=options.funcs.kind;
else
    name=pts(1).kind;
end
if isempty(options.lgname)
    lgname=name;
else
    lgname=options.lgname;
end
if isempty(options.color)
    clarray=colormap('lines');
    cl=clarray(strcmp(name,typenames),:);
else
    cl=options.color;
end
stfac=options.stability;
stcolor=1-(1-repmat(cl,length(nunst_unique),1)).*repmat(stfac.^nunst_unique(:),1,3);
%% divide up the branch into parts between special points
[bd,bifind,labels,markers]=splitbranch(pts,nunst,typenames,options.markersymbols);
nparts=size(bd,2);
nunstparts=nunst(floor(mean(bd,1)));
[nunstparts,ix]=sort(nunstparts,'descend');
bd=bd(:,ix);
%% plot bifurcation diagram
ax=options.ax;
ish=ishold(ax);
pdeco={'linewidth',options.linewidth,'markersize',options.markersize};
lg=options.oldlegend;
for i=1:nparts
    col=find(nunstparts(i)==nunst_unique,1,'first');
    lgline=plot(ax,p1(bd(1,i):bd(2,i)),p2(bd(1,i):bd(2,i)),'.-','color',stcolor(col,:),...
        pdeco{:});
    lgtext=stprint(lgname,nunstparts(i));
    set(lgline,'DisplayName',lgtext,'UserData','Plot2dBranch');
    lg{1}=[lg{1}(:)',lgline(:)'];
    lg{2}=[lg{2}(:)',{lgtext}];
    if i==1
        hold(ax,'on');
    end
end
%% plot codimension 2 bifurcation points
if ~isempty(markers)
    for i=1:length(bifind)
        lgpt=plot(ax,p1(bifind(i)),p2(bifind(i)),markers(i),...
            'markerfacecolor',options.bifcolor,'markeredgecolor',options.bifcolor,...
            pdeco{:});
        lgtext=sprintf('%s',labels{i});
        set(lgpt,'DisplayName',lgtext,'UserData','Plot2dBranch');
        lg{1}=[lg{1},lgpt(:)'];
        lg{2}=[lg{2},{lgtext}];
    end
end
%% make legend texts unique
obj=get(ax,'Children');
sel=arrayfun(@(c)~isempty(get(c,'UserData'))&&...
    strcmp(get(c,'UserData'),'Plot2dBranch'),obj);
obj=obj(sel);
lgtxt=get(obj,{'DisplayName'});
[lgnames,ix]=unique(lgtxt);
lghandles=obj(ix);
lg={lghandles,lgnames};
legend(ax,lg{1},lg{2});
if ~ish
    hold(ax,'off');
end
end
%%
function s=stprint_loc(name,i,maxi)
if maxi<=0
    s=sprintf('%s',name);
elseif i<maxi
    s=sprintf('%s #unst=%d',name,i);
else
    s=sprintf('%s #unst>=%d',name,i);
end
end
%%
function [bd,bifind,labels,markers]=splitbranch(pts,nunst,types,symbols)
np=length(pts);
if isfield(pts(1),'flag') && ~isempty([pts.flag]) 
    % flags indicate computation of bifurcations occured
    bifind=find(arrayfun(@(x)~isempty(x.flag),pts));
    labels={pts(bifind).flag};
else
    bifind=find(diff(nunst(:)')~=0);
    labels={};
    markers={};
end
bd=unique([1,bifind,np]);
bd=[bd(1:end-1);bd(2:end)];
%%  assign markers to bifurcations
[biftypes,~,biflocations]=unique(labels);
[~,type_ind]=ismember(biftypes,types);
if ~isempty(labels)
    markers=symbols(mod(type_ind-1,length(symbols))+1);
    markers=markers(biflocations);
end
end
