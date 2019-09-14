function [branch,suc]=correct_ini(funcs,branch,pini,dir,step,correc)
%% set and correct first (two) initial points along branch
%
% used in SetupFold, SetupHopf, SetupPOFold, SetupTorusbifurcation,
% SetupPeriodDoubling
%
% [branch,suc]=correct_ini(funcs,branch,pini,dir,step,correc)
%
% inputs
% funcs:  (struct) functions structure defining problems and rhs
% branch: (stuct) branch structure with initially empty point field
% pini:   (struct) initial guess for first point on branch
% dir:    (int) index of parameter value to change for second point 
%         (this routine cannot be used for branching off),
% step:   (double) deviation taken in parameter dir for second point
% correc: (logical) whether to perform Newton correction or not
%
% $Id: correct_ini.m 369 2019-08-27 00:07:02Z jansieber $
%
%%
suc=true;
mth=branch.method.point;
if correc
    %% correct initial guess
    free_par=parameter_fix(branch.parameter,pini);
    pnull=p_tangent(funcs,mth,pini,free_par);
    [pfirst,suc]=p_correc(funcs,pini,free_par,...
        pnull,mth,0,pini);
    if ~suc
        pfirst=pini;
        warning('correct_ini:fail','Correction failed');
    end
else
    pfirst=pini;
end
branch=rmfield(branch,'point');
branch.point(1)=pfirst;
if ~suc
    return
end
%% append 2nd point if desired
if ~isempty(dir)
    p2=branch.point(1);
    free_par=branch.parameter.free;
    pnull=p_tangent(funcs,mth,p2,free_par);
    assert(~isempty(pnull),'correct_ini: empty tangent space');
    itan=1;
    psgn=1;
    if dir(1)~=0
        psgn=arrayfun(@(p)p.parameter(dir(1)),pnull);
        [mx,itan]=max(abs(psgn));
        if mx==0
            warning('correct_ini:step',...
                'tangent has no component in parameter direction %d',dir(1));
            itan=1;
            psgn=1;
        end
    end
    p2=axpy(step*sign(psgn(itan)),pnull(itan),p2,free_par);
    if correc
        [p2c,suc]=p_correc(funcs,p2,free_par,pnull,mth,0,p2);
        if suc
            branch.point(2)=p2c;
        else
            warning('correct_ini:fail','Correction failed');
        end
    else
        p2c=p2;
    end
    branch.point(2)=p2c;
end
end
%% Determine if free parameter is on boundary and remove from set of free parameters
function free_par=parameter_fix(bpar,point)
free_par=bpar.free;
val=point.parameter(free_par);
bounds=cat(1,bpar.max_bound,bpar.min_bound);
for i=1:size(bounds,1)
    ind=find(bounds(i,1)==free_par);
    if isempty(ind)
        continue
    end
    if val(ind)==bounds(i,2)
        free_par(ind)=[];
    end
    if isempty(free_par)
        break
    end
end
end
%% add two points
function pnew=axpy(a,p1,p2,free_par)
x1=dde_x_from_point(p1,free_par);
x2=dde_x_from_point(p2,free_par);
pnew=dde_point_from_x(a*x1+x2,p1,free_par);
end
