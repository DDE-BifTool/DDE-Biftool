function p=p_axpy(a,x,y)

% function p=p_axpy(a,x,y)
% INPUT:
%	a (real) scalar 
%	x first point
%	y second point 
% OUTPUT:
%	p resulting point a*x+y

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

if isempty(a) 
  error('P_AXPY: empty scalar!');
elseif isempty(x) 
  error('P_AXPY: empty x!');
end;

p.kind=x.kind;

switch x.kind,
  case 'hcli',
    if isempty(y)
      p.parameter=a*x.parameter;
      p.mesh=x.mesh;
      p.degree=x.degree;
      p.profile=a*x.profile;
      p.period=a*x.period;
      p.x1=a*x.x1;
      p.x2=a*x.x2;
      p.lambda_v=a*x.lambda_v;
      p.lambda_w=a*x.lambda_w;
      p.v=a*x.v;
      p.w=a*x.w;
      p.alpha=a*x.alpha;
      p.epsilon=x.epsilon;
    else
      p.parameter=a*x.parameter+y.parameter;
      p.mesh=x.mesh;
      p.degree=x.degree;
      if ~isempty(x.mesh)
        if isempty(y.mesh)
          ly=size(y.profile,2);
          p.profile=a*x.profile+hcli_eva(y.profile,0:1/(ly-1):1,x.mesh,y.degree);
        else
          p.profile=a*x.profile+hcli_eva(y.profile,y.mesh,x.mesh,y.degree);
        end;
      else
        lx=size(x.profile,2);
        if ~isempty(y.mesh) 
          p.profile=a*x.profile+hcli_eva(y.profile,y.mesh,0:1/(lx-1):1,y.degree);
        else
          ly=size(y.profile,2);
          if lx~=ly
            p.profile=a*x.profile+hcli_eva(y.profile,0:1/(ly-1):1,0:1/(lx-1):1,y.degree);
          else
            p.profile=a*x.profile+y.profile; 
          end;
        end;
      end;
      p.period=a*x.period+y.period;
      p.x1=a*x.x1+y.x1;
      p.x2=a*x.x2+y.x2;
      p.lambda_v=a*x.lambda_v+y.lambda_v;
      p.lambda_w=a*x.lambda_w+y.lambda_w;
      p.v=a*x.v+y.v;
      p.w=a*x.w+y.w;
      p.alpha=a*x.alpha+y.alpha;
      % extra line to normalize alpha's
      if norm(p.alpha)~=0,
        p.alpha=p.alpha/norm(p.alpha);
      end;
      p.epsilon=x.epsilon;
    end;

  case 'psol',
    if isempty(y)
      p.parameter=a*x.parameter;
      p.mesh=x.mesh;
      p.degree=x.degree;
      p.profile=a*x.profile;
      p.period=a*x.period;
    else
      p.parameter=a*x.parameter+y.parameter;
      p.mesh=x.mesh;
      p.degree=x.degree;
      if ~isempty(x.mesh)
        if isempty(y.mesh)
          ly=size(y.profile,2);
          p.profile=a*x.profile+psol_eva(y,x.mesh);
        else
          p.profile=a*x.profile+psol_eva(y,x.mesh);
        end;
      else
        lx=size(x.profile,2);
        if ~isempty(y.mesh) 
          p.profile=a*x.profile+psol_eva(y.profile,y.mesh,0:1/(lx-1):1,y.degree);
        else
          ly=size(y.profile,2);
          if lx~=ly
            p.profile=a*x.profile+psol_eva(y.profile,0:1/(ly-1):1,0:1/(lx-1):1,y.degree);
          else
            p.profile=a*x.profile+y.profile; 
          end;
        end;
      end;
      p.period=a*x.period+y.period;
    end;
  case 'hopf',
    if isempty(y)
      p.parameter=a*x.parameter;
      p.x=a*x.x;
      p.v=a*x.v;
      p.omega=a*x.omega;
    else
      p.parameter=a*x.parameter+y.parameter;
      p.x=a*x.x+y.x;
      p.v=a*x.v+y.v;
      p.omega=a*x.omega+y.omega;
    end;
  case 'fold',
    if isempty(y)
      p.parameter=a*x.parameter;
      p.x=a*x.x;
      p.v=a*x.v;
    else
      p.parameter=a*x.parameter+y.parameter;
      p.x=a*x.x+y.x;
      p.v=a*x.v+y.v;
    end;
  case 'stst',
    if isempty(y)
      p.parameter=a*x.parameter;
      p.x=a*x.x;
    else
      p.parameter=a*x.parameter+y.parameter;
      p.x=a*x.x+y.x;
    end;
  otherwise, 
    err=x.kind,
    error('P_AXPY: point is not recognized!');
end;
%% all other fields present in x are initialized in p as empty
extrafields=setdiff(fieldnames(x),fieldnames(p));
for i=1:length(extrafields)
    p.(extrafields{i})=[];
end
p=orderfields(p,x);

end
