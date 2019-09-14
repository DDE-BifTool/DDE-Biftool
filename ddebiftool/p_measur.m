function s=p_measur(p,m)

% function sc_measure=p_measur(p,measure)
% INPUT:
%	p point
%	measure measure struct
% OUTPUT:
%	sc_measure scalar measure of point

% (c) DDE-BIFTOOL v. 1.00, 06/05/2000

row=m.row;
col=m.col;

if length(row)==3
    row(4)=' ';
end
if length(col)==3
    col(4)=' ';
end

f=getfield(p,m.field); %#ok<*GFLD>

if isempty(f)
    s=[];
    return;
end

if ~isempty(m.subfield)
    f=getfield(f,m.subfield);
    if isempty(f)
        s=[];
        return;
    end
end

switch col
    case {'max','max '}
        s=max(f(row,:));
    case {'min','min '}
        s=min(f(row,:));
    case 'mean'
        s=mean(f(row,:));
    case 'ampl'
        s=max(f(row,:))-min(f(row,:));
    otherwise
        switch row
            case {'max ','max'}
                s=max(f(:,col));
            case {'min ','min'}
                s=min(f(:,col));
            case 'mean'
                s=mean(f(:,col));
            case 'ampl'
                s=max(f(:,col))-min(f(:,col));
            case {'all','all '}
                switch col
                    case {'all','all '}
                        s=f(:,:);
                    otherwise
                        s=f(:,col);
                end
            otherwise
                switch col
                    case {'all','all '}
                        s=f(row,:);
                    otherwise
                        s=f(row,col);
                end
        end
end
if ~isempty(m.func)
  s=feval(m.func,s);
end
end
