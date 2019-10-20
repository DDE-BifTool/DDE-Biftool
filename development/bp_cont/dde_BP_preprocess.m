function data=dde_BP_preprocess(data)
if ~isfield(data,'pass_on')
    data.pass_on={};
end
f=@(x)p_correc_rhs(data.funcs,data.method,data.point,data.free_par,...
    'x',x,'pref',data.previous,'step_cond',data.step_cnd,...
    'extracolumns',data.extracolumns,data.pass_on{:});
x0=dde_x_from_point(data.point,data.free_par);
[r,J]=f(x0); %#ok<ASGLU>
data.nextrapar=data.funcs.ip.nextrapar;
if data.nextrapar>0
    data.extracolumns=dde_nullspaces_lr(J','nulldim',data.nextrapar);
    data.nuserpar=length(data.point.parameter);
    %data.point.parameter=cat(2,data.point.parameter,zeros(1,data.nextrapar));
    %data.previous.parameter=cat(2,data.previous.parameter,zeros(1,data.nextrapar));
    %data.free_par=cat(1,data.free_par,data.nuserpar+(1:data.nextrapar));
end
end
    