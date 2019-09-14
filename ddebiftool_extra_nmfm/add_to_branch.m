function [newbranch,success] = add_to_branch(funcs,newbranch,branch,i,bif_point)

new_point=branch.point(i);
switch new_point.kind
    case 'stst'
        new_degpoint=dde_stst_create('point',bif_point);
    case 'fold'
        new_degpoint=p_tofold(bif_point,'funcs',funcs);
    case 'hopf'
        new_degpoint=p_tohopf(bif_point,'funcs',funcs);
end

rtolfactor = branch.method.bifurcation.radial_tolerance_factor;

old_point = branch.point(i-1);
oldtonew = p_axpy(-1, old_point, new_point); % new_point - old_point
direction = p_axpy(1/p_norm(oldtonew),oldtonew,[]); % normalize
oldtodeg = p_axpy(-1, old_point, new_degpoint); % new_degpoint - old_point
projection = p_axpy(p_inprod(oldtodeg, direction),direction,old_point);
oldtoproj = p_axpy(-1,old_point,projection);
newtoproj = p_axpy(-1, oldtoproj, oldtodeg);
% Hack to disregard v
newtoproj.v = 0;
radial_distance = p_norm(newtoproj);
radiustol = rtolfactor*p_norm(oldtonew);

halfway_point = p_axpy(1,new_point,branch.point(i-1));
halfway_point = p_axpy(1/2,halfway_point,[]);

line_distance = p_norm(p_axpy(-1,halfway_point, projection));
linetol = p_norm(oldtonew)/2;

if ~(line_distance <= linetol && radial_distance <= radiustol)
    fprintf('BR_BIFDET: the detected %s point does not fall within the branch.\n\n', bif_point.kind);
    success = 0;
else
    if isfield(bif_point,'stability')
        new_degpoint.stability=bif_point.stability;
    end
    if ~isfield(newbranch.point(end),'stability') && isfield(new_degpoint,'stability')
        new_degpoint=rmfield(new_degpoint,'stability');
    end
    new_degpoint.flag=bif_point.kind;
    if isfield(bif_point,'nmfm')
        new_degpoint.nmfm=bif_point.nmfm;
    end
    if strcmp(new_degpoint.kind,'BT')
        new_degpoint = rmfield(new_degpoint,'p1');
        new_degpoint = rmfield(new_degpoint,'p0');
    end
    newbranch.point = [newbranch.point(1:i-1), new_degpoint, newbranch.point(i:end)];
    success=1;
end

end