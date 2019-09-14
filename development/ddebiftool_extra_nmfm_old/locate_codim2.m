function [new_bifpoint,success] = locate_codim2(funcs,branch,point,ind,bifkind)
max_iter = branch.method.bifurcation.secant_iterations;
free_par = branch.parameter.free;
method = branch.method.point;
stmethod = branch.method.stability;
conv_r = branch.method.bifurcation.secant_tolerance;
imagthresh=branch.method.bifurcation.imagthreshold; % threshold for treating a root as complex
isimag=@(x)x>imagthresh;
isreal=@(x)~isimag(x);

switch bifkind
    case 'cusp'
        cursign = point.nmfm.b; 
    case 'ZH'
        cursign = nmfm_smrp(funcs, point, stmethod,...
            'remove_omega',true,'threshold',isreal);
    case 'ZH_f'
        cursign = nmfm_smrp(funcs, point, stmethod,...
            'remove_omega',false,'threshold',isimag);
    case 'HH'
        cursign = nmfm_smrp(funcs, point, stmethod,...
            'remove_omega',true,'threshold',isimag);
    case 'GH'
        cursign = point.nmfm.L1;
    case 'BT' % use definig system
        [new_bifpoint,success]=locate_BT(funcs,branch,ind);
        return;
end

prevpoint=branch.point(ind-1);
if cursign > 0
    pospoint = point;
    negpoint = prevpoint;
else
    pospoint = prevpoint;
    negpoint = point;
end
secant=p_axpy(-1,point,prevpoint);
secant=p_secant(secant,p_norm(point));
success=0;
for ind = 1:max_iter
    sumpoint = p_axpy(1, pospoint, negpoint);
    halfpoint = p_axpy(0.5, sumpoint, []);
    halfpoint = p_correc(funcs, halfpoint,free_par,secant,method);
    switch bifkind
        case 'cusp'
            halfpoint = nmfm_fold(funcs,halfpoint);
            cursign = halfpoint.nmfm.b; 
        case 'ZH'
            [cursign, stability] = nmfm_smrp(funcs, halfpoint, stmethod,...
                'remove_omega',true,'threshold',isreal);
        case 'ZH_f'
            [cursign, stability] = nmfm_smrp(funcs, halfpoint, stmethod,...
                'remove_omega',false,'threshold',isimag);
        case 'HH'
            [cursign, stability] = nmfm_smrp(funcs, halfpoint, stmethod,...
                'remove_omega',true,'threshold',isimag);
        case 'GH'
            halfpoint = nmfm_hopf(funcs,halfpoint, prevpoint);
            cursign = halfpoint.nmfm.L1;
    end
    if abs(cursign) < conv_r
        if exist('stability')
            halfpoint.stability=stability;
        end
        success = 1;
        break;
    end
    if cursign < 0
        negpoint = halfpoint;
    else
        pospoint = halfpoint;
    end
end
if success
    switch bifkind
        case 'cusp'
            new_bifpoint = p_tocusp(halfpoint);
        case 'ZH'
            new_bifpoint = p_tozeho(halfpoint);
        case 'ZH_f'
            new_bifpoint = p_tozeho(halfpoint);
        case 'HH'
            new_bifpoint = p_tohoho(halfpoint);
        case 'GH'
            new_bifpoint = p_togenh(halfpoint);
    end
else
    new_bifpoint=[];
end
end
