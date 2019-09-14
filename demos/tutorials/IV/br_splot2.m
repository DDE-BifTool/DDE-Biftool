function br_splot2(br,nunst,par,ordinate,varargin)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: br_splot2.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
    pars=arrayfun(@(x)x.parameter(par),br.point);
    clrs=colormap('lines');
    pernunst_cases=[0; find(diff(nunst)~=0); length(nunst)-1];
    
    if strcmp(br.point(1).kind,'psol')
        pprofs1=cell2mat(arrayfun(@(x)x.profile(1,:)',br.point,...
            'uniformoutput',false));
        if strcmp(ordinate,'amplitude')
            amp=max(pprofs1)-min(pprofs1);
            ylabel_str='amplitude';
        else
            amp=max(pprofs1);
            ylabel_str='$\max x_1$';
        end
    else
        amp=arrayfun(@(p)p.x(1),br.point);
        ylabel_str='$x_1$';
    end

    xlabel('$c$','Interpreter','LaTex');
    ylabel(ylabel_str,'Interpreter','LaTex');

    % set line width if requested
    if nargin==5
        lw=varargin{1};
    else
        lw=1.2;
    end
    
    if strcmp(br.point(1).kind,'stst')
        bgetpar=@(x,i,bif)arrayfun(@(p)p.parameter(i),x.point(br_getflags(x,bif)));
        bgetx=@(x,i,bif)arrayfun(@(p)p.x(i),x.point(br_getflags(x,bif)));
        % strange fix for GNU Octave
        % without the first plot an error is thrown
        % need to look into this
        hb1=plot(bgetpar(br,par,'hopf'),bgetx(br,1,'hopf'),'ks','MarkerSize',4);
        hold on
        hb1=plot(bgetpar(br,par,'hopf'),bgetx(br,1,'hopf'),'ks','MarkerSize',4);
        hb2=plot(bgetpar(br,par,'fold'),bgetx(br,1,'fold'),'ro','MarkerSize',4);
    else  % add dots at bifurcation points          
        for i=2:length(pernunst_cases)-1
            plot(pars(pernunst_cases(i)+1),amp(pernunst_cases(i)+1),...
                '.','color','k');
        end
        hold on
    end

    for i=1:length(pernunst_cases)-1
        sel=pernunst_cases(i)+1:pernunst_cases(i+1)+1;
        h(nunst(pernunst_cases(i+1))+1)=plot(pars(sel),amp(sel),...
            '-','color',clrs(nunst(pernunst_cases(i+1))+1,:),...
            'LineWidth',lw);
    end
    
    % add legend
    nunst_unique=unique(nunst);
    lstr={};
    for i=1:length(nunst_unique)
        lstr={lstr{:},sprintf('#unst=%d',nunst_unique(i))};
    end
    
    
    if strcmp(br.point(1).kind,'psol')
        legend(h(nunst_unique+1),lstr{:})
    else
        if isempty(hb2)
            legend([h(nunst_unique+1) hb1],{lstr{:},'Hopf'})
        else
            legend([h(nunst_unique+1) hb1 hb2],{lstr{:},'Hopf','fold'})
        end
    end
end
