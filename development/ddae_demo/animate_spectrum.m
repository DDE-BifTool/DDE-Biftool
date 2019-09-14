function ax=animate_spectrum(branch,ip,varargin)
default={'figures',[3,4,5],'order','','pause',Inf,...
    'xlim',[-3,2],'ylim',8*[-1,1],'hold','off','ip_essential',[],...
    'nunst',[],'efcomp',1,'solcomp',1};
options=dde_set_options(default,varargin,'pass_on');
lw={'linewidth',2};
pt=branch.point(1);
if isfield(pt,'x')
    field='l0';
    isfloquet=false;
else
    field='mu';
    isfloquet=true;
end
[ipvals,ipnames]=parnames(branch.parameter.free,ip);
getpar=@(i)arrayfun(@(p)p.parameter(i),branch.point);
xbif=getpar(ipvals(1));
ybif=getpar(ipvals(2));
switch options.order
    case 'reverse'
        seq=length(branch.point):-1:1;
    otherwise
        seq=1:length(branch.point);
end
figure(options.figures(1));clf;ax1=gca;hold(ax1,options.hold);
if length(options.figures)>=2
    plotbif=true;
    figure(options.figures(2));clf;ax2=gca;
    plot(ax2,xbif,ybif,'.-',lw{:});
    hold(ax2,'on');
else
    plotbif=false;
end
th=linspace(0,2*pi,100);
for i=seq
    pt=branch.point(i);
    evl=pt.stability.(field);
    plot(ax1,real(evl),imag(evl),'o','linewidth',2);
    grid(ax1,'on');
    if ~isfloquet
        ax1.XLim=options.xlim;
        ax1.YLim=options.ylim;
    else
        hold(ax1,'on');
        if ~isempty(options.ip_essential)
            rad=pt.parameter(ip.(options.ip_essential));
            plot(ax1,rad*cos(th),rad*sin(th),lw{:});
            nmu=find(abs(evl)>=rad);
        else
            nmu=length(evl);
        end
        plot(ax1,real(evl),imag(evl),'o',cos(th),sin(th),lw{:});
        axis(ax1,'equal');
        if length(options.figures)>=3
            figure(options.figures(3));
            subplot(2,1,1);ax3=gca;
            subplot(2,1,2);ax4=gca;
            ef=pt.stability.eigenfuncs;
            hold(ax3,'off');
            for k=1:length(nmu)
                plot(ax3,ef(nmu(k)).mesh*pt.period,...
                real(ef(nmu(k)).profile(options.efcomp,:)),lw{:});
                hold(ax3,'on');
            end
            hold(ax3,'off');
            title(ax3,sprintf('Eigenfunctions comp %d',options.efcomp))
            plot(ax4,pt.mesh*pt.period,pt.profile(options.solcomp,:),lw{:});
            title(ax4,sprintf('Solution comp %d',options.solcomp))
        end
    end
    stc=cellfun(@(a,b)sprintf('%s=%g ',a,pt.parameter(b)),ipnames,...
        num2cell(ipvals),'uniformoutput',false);
    if ~isempty(options.nunst)
        stc{end+1}=sprintf(', nunst=%d',options.nunst(i));
    end
    title(ax1,[stc{:}]);
    hold(ax1,options.hold);
    if plotbif && i==seq(1)
        xybif=plot(ax2,xbif(i),ybif(i),'o',lw{:});
    else
        xybif.XData=xbif(i);
        xybif.YData=ybif(i);
    end
    drawnow
    locpause(options.pause);
end
end
%%
function [ind,names]=parnames(parfree,ip)
ipnames=fieldnames(ip);
ipvals=cell2mat(struct2cell(ip));
ind=intersect(parfree,ipvals);
names=ipnames(ind);
end
%%
function locpause(ptime)
if isinf(ptime)
    pause
elseif ptime>0
    pause(ptime);
end
end

