function pt=dde_stabilitykind(pt,data,args)
%% determine and reduce point type, determine stability
% $Id: dde_stabilitykind.m 343 2019-05-12 00:01:02Z jansieber $
%%
if ~isstruct(pt)
    pt=dde_point_from_x(pt,data.info.point,data.ipar);
end
kind=pt.kind;
funcs=data.funcs;
if strcmp(kind,'psol') && isfield(funcs,'kind')
    kind=funcs.kind;
    pt=funcs.get_comp(pt,'solution_for_stability');
    funcs=funcs.userfuncs;
end
if iscell(args)
    pt.stability=p_stabil(funcs,pt,data.info.method.stability,args{:});
end
pt.kind=kind;
end
