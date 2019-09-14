function [branch,succ,fail]=br_recmp(funcs,branch,point_numbers)
%% recompute branch partially at selected points
% function [recmp_branch,succ,fail]=br_recmp(funcs,branch,point_numbers)
% INPUT:
%   funcs problem functions
%	branch 
%	point_numbers numbers of points to recompute or [] for all
% OUTPUT:
%	recmp_branch (partually or completely) recomputed branch
%	succ number of successfull corrections
%	fail number of failed corrections

% (c) DDE-BIFTOOL v. 2.00, 23/12/2000
%
% $Id: br_recmp.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
free_par=branch.parameter.free;

ll=length(branch.point);

if ll<1 
    error('BR_RECMP: branch is empty: ll=%d points!',ll);
end;

if ~exist('point_numbers','var')
    point_numbers=[];
end;

if isempty(point_numbers)
    point_numbers=1:ll;
end;

succ=0;

for i=1:length(point_numbers) 
    j=point_numbers(i);
    if j>1 
        left=branch.point(j-1);
    else
        left=branch.point(j);
    end;
    if j<ll
        right=branch.point(j+1);
    else
        right=branch.point(j);
    end;
    secant=p_axpy(-1,right,left);  

    if strcmp(secant.kind,'psol') || strcmp(secant.kind,'hcli')
        if length(secant.mesh)==length(branch.point(j).mesh) && secant.degree==branch.point(j).degree,
            secant=p_secant(secant,p_norm(branch.point(j)));
            [p,success]=p_correc(funcs,branch.point(j),free_par,secant,branch.method.point,j+1);
        else
         if isempty(branch.point(j).mesh)
            msh=0:1/(size(branch.point(j).profile,2)-1):1;
         else
            msh=branch.point(j).mesh;
         end;
         secant=p_remesh(secant,branch.point(j).degree,msh);
         if isempty(branch.point(j).mesh)
             secant.mesh=[];
         end;
         secant=p_secant(secant,p_norm(branch.point(j)));
         [p,success]=p_correc(funcs,branch.point(j),free_par,secant,branch.method.point,j+1);    
        end;
    else
        [p,success]=p_correc(funcs,branch.point(j),free_par,secant,branch.method.point,j+1); 
    end;

    if success  
        succ=succ+1;
        branch.point(j)=p;    
    else      
        s=sprintf('BR_RECMP warning: failure during recomputation of point %d.',j);
        disp(s);
    end;
end;

fail=length(point_numbers)-succ;

return;
