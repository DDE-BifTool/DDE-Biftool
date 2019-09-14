function out=dde_psolfold_extract(point,comp)
%% extract components of psolfold
% $Id: dde_psolfold_extract.m 338 2019-05-09 17:35:35Z jansieber $
%%
if nargin<1
    out={'solution','var','delays'};
    return
end
if isempty(point)
    out=repmat(dde_psol_create(),size(point));
    return
end
switch comp
    case {'solution','psol'}
        out=arrayfun(@(p)dde_psol_create('point',p),point);
    case {'var','variation'}
        out=arrayfun(@(p)dde_psol_create('point',p),point);
        for i=1:length(point)
            out(i).profile=point(i).profile_var;
            out(i).period=point(i).parameter_var(1)*point(i).period;
            out(i).parameter=0*point(i).parameter;
            out(i).parameter(point(i).free_par)=point(i).parameter_var(2:end);
        end
    case {'delays','delay'}
        out=arrayfun(@(p)p.parameter_delay(:),point,'uniformoutput',false);
        out=cell2mat(out);
end     
end