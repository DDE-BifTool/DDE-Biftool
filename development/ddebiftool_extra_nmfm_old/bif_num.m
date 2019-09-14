function result = bif_num(pointtype,tofrom)
%% Convert bifurcation type to its index
% Return number of supported types on 'count'
%
% $Id: bif_num.m 114 2015-09-02 02:24:46Z jansieber $
%%
typelist={...
    '',...     % no special point
    'stst',... % steady state
    'hopf',... % Hopf bifurcation
    'fold',... % fold bifurcation
    'psol',... % periodic orbit
    'hcli',... % connecting orbit
    'genh',... % generalized Hopf
    'hoho',... % double Hopf (Hopf-Hopf)
    'zeho',... % zero-Hopf (Gavrilov-Guckenheimer)
    'BT',...   % Takens-Bogdanov bifurcation
    'CP',...   % cusp
    'count'};
if nargin<2
%     fprintf('known point types:\n');
%     for i=1:length(typelist)-1
%         fprintf('type %2d: ''%s''\n',i-2,typelist{i});
%     end
    result=typelist(1:end-1);
    return
end
switch tofrom
    case {'->','to','To'}
        if ischar(pointtype)
            pointtype={pointtype};
        end
        result=NaN(size(pointtype));
        for i=1:numel(result)
            num=find(strcmp(pointtype{i},typelist));
            if isempty(num)
                warning('bif_num:unknown','type %s is unknown',pointtype{i});
            else
                result(i)=num-2;
            end
        end
        % adjust count
        result(result==length(typelist)-2)=result(result==length(typelist)-2)-1;
    case {'<-','from','From'}
        result=cell(size(pointtype));
        for i=1:numel(result)
            if pointtype(i)>=-1 && pointtype(i)<length(typelist)-2
                result(i)=typelist(pointtype(i)+2);
            else
                warning('bif_num:unknown','type %d is unknown',pointtype(i));                
            end
        end
        if numel(result)==1
            result=result{1};
        end
    otherwise
        error('bif_num:arg','illegal 2nd argument %s',tofrom);
end
end
