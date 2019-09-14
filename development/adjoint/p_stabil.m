function stability=p_stabil(funcs,p,method)
%% compute stability information for point
% function stability=p_stabil(funcs,point,method)
% INPUT:
%   funcs problem functions
%	point solution point
%	method method parameters 
% OUTPUT:
%	stability stability information

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
% Update on 05/03/2007 ("flag_newhheur"  <=> (imag(method.lms_parameter_rho)~=0) )   
%
% $Id: p_stabil.m 19 2014-04-11 14:15:36Z jan.sieber $
%
%%
if strcmp(p.kind,'psol')    
  if isempty(p.mesh)
    mesh=0:1/(size(p.profile,2)-1):1;
  else
    mesh=p.mesh;
  end;
  rho=method.minimal_modulus;
  max_number=method.max_number_of_eigenvalues;
  if isempty(method.collocation_parameters)
    col=poly_gau(p.degree);
  else
    col=method.collocation_parameters;
  end;
  if ~funcs.tp_del
    mu=mult_app(funcs,p,rho,max_number,col);
  else 
    d_ac=method.delay_accuracy;
    mu=mult_app(funcs,p,rho,max_number,col,d_ac);
  end;
  if ~isempty(mu)
    [ydum,index_vector]=sort(abs(mu)); %#ok<ASGLU>
    stability.mu=mu(index_vector(length(index_vector):-1:1));
  else
    stability.mu=[];
  end;
else
% p.kind=='stst' OR 'hopf' ...
  if imag(method.lms_parameter_rho)~=0,  
    stability=stst_stabil(funcs,p,method);
  else
    real_part=method.minimal_real_part;
    a=method.lms_parameter_alpha;
    b=method.lms_parameter_beta;
    ord=method.interpolation_order;
    h_min=method.minimal_time_step;
    h_max=method.maximal_time_step;
    rho=method.lms_parameter_rho;
    if ~funcs.tp_del
      [l0,h]=root_app(funcs,p.x,real_part,a,b,ord,p.parameter,h_min,h_max,rho);
    else
      d_ac=method.delay_accuracy;
      [l0,h]=root_app(funcs,p.x,real_part,a,b,ord,p.parameter,h_min,h_max,rho,d_ac);
    end;
    if length(l0)>method.max_number_of_eigenvalues
      l0=l0(1:method.max_number_of_eigenvalues);
    end;
    stability.h=h;
    stability.l0=l0;
    max_n=method.max_newton_iterations;
    epsi=method.root_accuracy;
    if max_n==0 || isempty(l0)
      stability.l1=[];
      stability.n1=[];
    else
      [l1,n1]=root_nwt(funcs,p.x,l0,max_n,epsi,p.parameter);
      if method.remove_unconverged_roots % remove n1 and sort roots
	l1=l1(n1~=-1);
	n1=[];
	if ~isempty(l1)
        [ydum,index_vector]=sort(real(l1)); %#ok<ASGLU>
        l1=l1(index_vector(length(index_vector):-1:1))';
	end;
      end;
      stability.l1=l1;
      stability.n1=n1;
    end;
  end;  
end;

return;
