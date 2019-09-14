function [data, y,J] = dde_rhs(~, data, u)
%% main r.h.s. for dde-biftool points
% $Id: dde_rhs.m 345 2019-05-12 00:11:15Z jansieber $
%%
f=p_correc_setup(data.funcs,data.info.point,data.ipar,data.info.method.point,...
    'remesh_flag',false,'previous',data.info.point);
[y,J]=f(u);
end