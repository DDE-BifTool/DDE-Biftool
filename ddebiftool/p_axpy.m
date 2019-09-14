function p=p_axpy(a,p_x,p_y,free_par)

% function p=p_axpy(a,x,y,free_par)
% INPUT:
%	a (real) scalar 
%	p_x first point
%	p_y second point 
% OUTPUT:
%	p resulting point a*p_x+p_y

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

if isempty(a) 
    error('P_AXPY: empty scalar!');
elseif isempty(p_x) 
    error('P_AXPY: empty x!');
end
if nargin<4
    free_par=1:length(p_x.parameter);
end
x=dde_x_from_point(p_x,free_par);
y=0*x;
if nargin>2 && ~isempty(p_y)
    p_y=p_remesh(p_y,p_x);
    y=dde_x_from_point(p_y,free_par);
end
axpy=a*x+y;
p=dde_point_from_x(axpy,p_x,free_par);
end
% p.kind=p_x.kind;
% 
% switch p_x.kind,
%   case 'hcli',
%     if isempty(p_y)
%       p.parameter=a*p_x.parameter;
%       p.mesh=p_x.mesh;
%       p.degree=p_x.degree;
%       p.profile=a*p_x.profile;
%       p.period=a*p_x.period;
%       p.x1=a*p_x.x1;
%       p.x2=a*p_x.x2;
%       p.lambda_v=a*p_x.lambda_v;
%       p.lambda_w=a*p_x.lambda_w;
%       p.v=a*p_x.v;
%       p.w=a*p_x.w;
%       p.alpha=a*p_x.alpha;
%       p.epsilon=p_x.epsilon;
%     else
%       p.parameter=a*p_x.parameter+p_y.parameter;
%       p.mesh=p_x.mesh;
%       p.degree=p_x.degree;
%       if ~isempty(p_x.mesh)
%         if isempty(p_y.mesh)
%           ly=size(p_y.profile,2);
%           p.profile=a*p_x.profile+hcli_eva(p_y.profile,0:1/(ly-1):1,p_x.mesh,p_y.degree);
%         else
%           p.profile=a*p_x.profile+hcli_eva(p_y.profile,p_y.mesh,p_x.mesh,p_y.degree);
%         end;
%       else
%         lx=size(p_x.profile,2);
%         if ~isempty(p_y.mesh) 
%           p.profile=a*p_x.profile+hcli_eva(p_y.profile,p_y.mesh,0:1/(lx-1):1,p_y.degree);
%         else
%           ly=size(p_y.profile,2);
%           if lx~=ly
%             p.profile=a*p_x.profile+hcli_eva(p_y.profile,0:1/(ly-1):1,0:1/(lx-1):1,p_y.degree);
%           else
%             p.profile=a*p_x.profile+p_y.profile; 
%           end;
%         end;
%       end;
%       p.period=a*p_x.period+p_y.period;
%       p.x1=a*p_x.x1+p_y.x1;
%       p.x2=a*p_x.x2+p_y.x2;
%       p.lambda_v=a*p_x.lambda_v+p_y.lambda_v;
%       p.lambda_w=a*p_x.lambda_w+p_y.lambda_w;
%       p.v=a*p_x.v+p_y.v;
%       p.w=a*p_x.w+p_y.w;
%       p.alpha=a*p_x.alpha+p_y.alpha;
%       % extra line to normalize alpha's
%       if norm(p.alpha)~=0,
%         p.alpha=p.alpha/norm(p.alpha);
%       end;
%       p.epsilon=p_x.epsilon;
%     end;
% 
%   case 'psol',
%     if isempty(p_y)
%       p.parameter=a*p_x.parameter;
%       p.mesh=p_x.mesh;
%       p.degree=p_x.degree;
%       p.profile=a*p_x.profile;
%       p.period=a*p_x.period;
%     else
%       p.parameter=a*p_x.parameter+p_y.parameter;
%       p.mesh=p_x.mesh;
%       p.degree=p_x.degree;
%       if ~isempty(p_x.mesh)
%         if isempty(p_y.mesh)
%           ly=size(p_y.profile,2);
%           p.profile=a*p_x.profile+psol_eva(p_y.profile,0:1/(ly-1):1,p_x.mesh,p_y.degree);
%         else
%           p.profile=a*p_x.profile+psol_eva(p_y.profile,p_y.mesh,p_x.mesh,p_y.degree);
%         end;
%       else
%         lx=size(p_x.profile,2);
%         if ~isempty(p_y.mesh) 
%           p.profile=a*p_x.profile+psol_eva(p_y.profile,p_y.mesh,0:1/(lx-1):1,p_y.degree);
%         else
%           ly=size(p_y.profile,2);
%           if lx~=ly
%             p.profile=a*p_x.profile+psol_eva(p_y.profile,0:1/(ly-1):1,0:1/(lx-1):1,p_y.degree);
%           else
%             p.profile=a*p_x.profile+p_y.profile; 
%           end;
%         end;
%       end;
%       p.period=a*p_x.period+p_y.period;
%     end;
%   case 'hopf',
%     if isempty(p_y)
%       p.parameter=a*p_x.parameter;
%       p.x=a*p_x.x;
%       p.v=a*p_x.v;
%       p.omega=a*p_x.omega;
%     else
%       p.parameter=a*p_x.parameter+p_y.parameter;
%       p.x=a*p_x.x+p_y.x;
%       p.v=a*p_x.v+p_y.v;
%       p.omega=a*p_x.omega+p_y.omega;
%     end;
%   case 'fold',
%     if isempty(p_y)
%       p.parameter=a*p_x.parameter;
%       p.x=a*p_x.x;
%       p.v=a*p_x.v;
%     else
%       p.parameter=a*p_x.parameter+p_y.parameter;
%       p.x=a*p_x.x+p_y.x;
%       p.v=a*p_x.v+p_y.v;
%     end;
%   case 'stst',
%     if isempty(p_y)
%       p.parameter=a*p_x.parameter;
%       p.x=a*p_x.x;
%     else
%       p.parameter=a*p_x.parameter+p_y.parameter;
%       p.x=a*p_x.x+p_y.x;
%     end;
%   otherwise, 
%     err=p_x.kind,
%     error('P_AXPY: point is not recognized!');
% end;
% %% all other fields present in x are initialized in p as empty
% extrafields=setdiff(fieldnames(p_x),fieldnames(p));
% for i=1:length(extrafields)
%     p.(extrafields{i})=[];
% end
% p=orderfields(p,p_x);
% 
% end
