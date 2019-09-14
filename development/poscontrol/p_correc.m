function [point,success]=p_correc(funcs,point,free_par,step_cnd,method,p_nr,...
    previous,d_nr,tz)
%% correct point using Newton iteration
% function [point,success]=p_correc(point0,free_par,step_cnd,method,adapt,
%                                 previous)
% INPUT:
%   funcs problem functions
%   point0 initial point guess
%   free_par free parameter numbers in N^d
%       step_cnd steplength condition(s) as point(s)
%   method method parameters
%   adapt if zero or absent, do not adapt mesh; if one, always adapt
%       previous (optional) previously computed branch point (used, in case
%            of periodic solutions or connecting orbits, to
%            minimize phase shift)
% OUTPUT:
%       point final solution
%       success nonzero for successfull computation
%
% OPTIONAL EXTRA INPUT:
%       d_nr number of delay crossing zero
%       tz (for psol only) time point for which tau(tz)=0 and dtau(tz)/dt=0
%
% (c) DDE-BIFTOOL v. 2.02, 16/6/2002
%
% $Id: p_correc.m 176 2017-03-13 00:25:33Z jansieber $
%
%%
sys_cond=funcs.sys_cond;


% optional parameters:
if nargin<=5
    p_nr=0;
end
if nargin<=6
    previous=point;
end
if nargin<=7
    d_nr=0;
end

% some method parameters:

max_iter=method.newton_max_iterations; % max number of newton iterations
nmon_iter=method.newton_nmon_iterations; % max number of nonmonotone iterations
conv_r=method.halting_accuracy; % required accuracy for convergence
print_r=method.print_residual_info; % print residual evolution

% replace old stability and other info by empty fields if present:
point=p_axpy(1,point,[]);

% initialize:

switch point.kind
    case 'stst'
        n=size(point.x,1);
        p_start=n;
    case 'fold'
        n=size(point.x,1);
        p_start=2*n;
        c=point.v'/(point.v'*point.v);
    case 'Hopf'
        n=size(point.x,1);
        p_start=3*n+1;
        vr=real(point.v);
        vi=imag(point.v);
        vn=vr'*vr+vi'*vi;
        c=point.v'/vn; % complex conjugate transpose
    case {'Psol','hcli'}
        n=size(point.profile,1);
        col=method.collocation_parameters;
        mm=point.degree;
        if ~isempty(col) && length(col)~=mm
            error('P_CORREC: number %d of collocation parameters differs from degree %d!',...
                mm,point.degree);
        end
        ll=floor((size(point.profile,2)-1)/mm);
        if ll~=(size(point.profile,2)-1)/mm
            error('P_CORREC: point mesh does not contain l intervals of m points!');
        end
        p_start=n*mm*ll+1+n;
        if d_nr==0
            ph=method.phase_condition;
        else
            ph=0;
        end;
        if isempty(point.mesh)
            mesh=0:1/(ll*mm):1;
        else
            mesh=point.mesh;
            ma=method.adapt_mesh_before_correct;
            if p_nr>1 && mod(p_nr,ma)==0
                % do not adapt mesh when p_nr=0
                % do not (yet) adapt mesh when p_nr=1
                % adapt mesh when p_nr>1 & mod(p_nr,ma)=0
                new_mesh=psol_msh(mesh,mm,point.profile,ll,point.degree);
                switch point.kind
                    case 'Psol'
                        point.profile=old_psol_eva(point.profile,mesh,new_mesh,mm);
                    case 'hcli'
                        point.profile=hcli_eva(point.profile,mesh,new_mesh,mm);
                    otherwise
                        pt_eva=str2func([point.kind,'_eva']);
                        point.profile=pt_eva(point,new_mesh);
                end
                point.mesh=new_mesh;
                mesh=new_mesh;
            end
        end
    otherwise
        %error('P_CORREC: point kind %s not recognized.',point.kind);
end

% store parameters:

par=point.parameter;

% perform Newton-Raphson iterations:

for i=1:max_iter
    
    % compute & apply corrections:
    resi=[];
    condi=repmat(point,0,1);
    if d_nr~=0
        [resi,condi]=delay_zero_cond(funcs,point,d_nr,free_par,tz);
    end
    if method.extra_condition
        [resi2,condi2]=sys_cond(point);
        resi=[resi;resi2(:)]; %#ok<AGROW>
        condi=[condi;condi2(:)]; %#ok<AGROW>
    end
    switch point.kind
        case 'stst'
            % produce jacobian
            [J,res]=stst_jac(funcs,point.x,par,free_par);
            % add (linear) steplength conditions
            for j=1:length(step_cnd)
                J(n+j,1:n)=step_cnd(j).x';
                for k=1:length(free_par)
                    J(n+j,n+k)=step_cnd(j).parameter(free_par(k));
                end;
                res(n+j,1)=0;
            end
            % add extra conditions
            res=[res;resi]; %#ok<AGROW>
            Jd=zeros(length(resi),size(J,2));
            for k=1:length(resi)
                Jd(k,:)=[condi(k).x',condi(k).parameter(free_par)];
            end
            J=[J;Jd]; %#ok<AGROW>
            % solve linear system
            if size(J,1)~=size(J,2)
                warning('p_correc:nonsquare','P_CORREC warning: use of nonsquare Jacobian.');
            end;
            
            dx=J\res;
            % apply non-parameter corrections
            point.x=point.x-dx(1:n);
            for j=1:length(free_par)
                par(free_par(j))=par(free_par(j))-real(dx(p_start+j));
            end;
            
            % fill in parameters:
            point.parameter=par;
            
        case 'fold'
            % produce jacobian
            [J,res]=fold_jac(funcs,point.x,point.v,par,free_par,c);
            % add (linear) steplength conditions
            for j=1:length(step_cnd)
                J(2*n+1+j,1:n)=step_cnd(j).x';
                J(2*n+1+j,n+1:2*n)=step_cnd(j).v';
                for k=1:length(free_par)
                    J(2*n+1+j,2*n+k)=step_cnd(j).parameter(free_par(k));
                end;
                res(2*n+1+j,1)=0;
            end;
            % add extra conditions
            res=[res;resi]; %#ok<AGROW>
            Jd=zeros(length(resi),size(J,2));
            for k=1:length(resi)
                Jd(k,:)=[condi(k).x',condi(k).v',condi(k).parameter(free_par)];
            end
            J=[J;Jd]; %#ok<AGROW>
            % solve linear system
            if size(J,1)~=size(J,2)
                warning('p_correc:nonsquare','P_CORREC warning: use of nonsquare Jacobian.');
            end;
            dx=J\res;
            % apply non-parameter corrections
            point.x=point.x-dx(1:n);
            point.v=point.v-dx(n+1:2*n);
            for j=1:length(free_par)
                par(free_par(j))=par(free_par(j))-real(dx(p_start+j));
            end;
            
            % fill in parameters:
            point.parameter=par;
           
        case 'Hopf'
            % produce jacobian
            [J,res]=hopf_jac(funcs,point.x,point.omega,point.v,par,free_par,c);
            % add (linear) steplength conditions
            for j=1:length(step_cnd)
                J(3*n+2+j,1:n)=step_cnd(j).x';
                J(3*n+2+j,n+1:2*n)=real(step_cnd(j).v)';
                J(3*n+2+j,2*n+1:3*n)=imag(step_cnd(j).v)';
                J(3*n+2+j,3*n+1)=step_cnd(j).omega;
                for k=1:length(free_par)
                    J(3*n+2+j,3*n+1+k)=step_cnd(j).parameter(free_par(k));
                end;
                res(3*n+2+j)=0;
            end;
            % add extra conditions
            res=[res;resi]; %#ok<AGROW>
            Jd=zeros(length(resi),size(J,2));
            for k=1:length(resi)
                Jd(k,:)=[condi(k).x',real(condi(k).v)',imag(condi(k).v)',...
                    condi(k).omega,condi(k).parameter(free_par)];
            end
            J=[J;Jd]; %#ok<AGROW>
            % solve linear system
            if size(J,1)~=size(J,2)
                warning('p_correc:nonsquare',...
                    'P_CORREC warning: use of nonsquare %dx%d Jacobian.',size(J,1),size(J,2));
            end;
            dx=J\res;
            % apply non-parameter corrections
            point.x=point.x-dx(1:n);
            point.v=point.v-dx(n+1:2*n)-sqrt(-1)*dx(2*n+1:3*n);
            point.omega=point.omega-dx(3*n+1);
            for j=1:length(free_par)
                par(free_par(j))=par(free_par(j))-real(dx(p_start+j));
            end;
            
            % fill in parameters:
            point.parameter=par;
            
        case 'Psol'
            % produce jacobian
            [J,res]=psol_jac(funcs,col,point.period,point.profile,mesh,mm,par,...
                free_par,ph);
            % add (linear) steplength conditions
            for j=1:length(step_cnd)
                for k=1:ll*mm+1
                    J(n*(mm*ll+1)+ph+j,(k-1)*n+1:k*n)=step_cnd(j).profile(:,k)';
                end;
                J(n*(mm*ll+1)+ph+j,n*(mm*ll+1)+1)=step_cnd(j).period;
                for k=1:length(free_par)
                    J(n*(mm*ll+1)+ph+j,n*(mm*ll+1)+1+k)=...
                        step_cnd(j).parameter(free_par(k));
                end;
                res(n*(mm*ll+1)+ph+j)=0;
            end;
            Jd=zeros(length(resi),size(J,2));
            for k=1:length(resi)
                Jd(k,:)=[condi(k).profile(:)',...
                    condi(k).period,condi(k).parameter(free_par)];
            end
            res=[res;resi]; %#ok<AGROW>
            J=[J;Jd]; %#ok<AGROW>
            % solve linear system
            if size(J,1)~=size(J,2)
                warning('p_correc:nonsquare','P_CORREC warning: use of nonsquare Jacobian.');
            end;
            dx=J\res;
            % apply non-parameter corrections
            for k=1:ll*mm+1
                point.profile(:,k)=point.profile(:,k)-dx((k-1)*n+1:k*n);
            end;
            point.period=point.period-dx(n*mm*ll+1+n);
            % apply parameter corrections:
            for j=1:length(free_par)
                par(free_par(j))=par(free_par(j))-real(dx(p_start+j));
            end;
            
            % fill in parameters:
            point.parameter=par;
            for j=1:length(free_par)
                par(free_par(j))=par(free_par(j))-real(dx(p_start+j));
            end;
            
            % fill in parameters:
            point.parameter=par;
            

        case 'hcli'
            % produce jacobian
            [J,res]=hcli_jac(funcs,col,point.period,point.profile,...
                mesh,mm,par,free_par,ph,...
                point.lambda_v,point.lambda_w,point.v,...
                point.w,point.alpha,point.epsilon,point.x1,...
                point.x2,previous);
            % add (linear) steplength conditions
            for j=1:length(step_cnd)
                sJ=size(J,1);
                for k=1:ll*mm+1
                    J(sJ+1,(k-1)*n+1:k*n)=step_cnd(j).profile(:,k)';
                end;
                J(sJ+1,n*(mm*ll+1)+1)=step_cnd(j).period;
                nb_par=length(free_par);
                for k=1:nb_par
                    J(sJ+1,n*(mm*ll+1)+1+k)=step_cnd(j).parameter(free_par(k));
                end;
                J(sJ+1,n*(mm*ll+1)+1+nb_par+(1:n))=step_cnd(j).x1';
                J(sJ+1,n*(mm*ll+1)+nb_par+n+1+(1:n))=step_cnd(j).x2';
                s1=length(step_cnd(j).v(1,:));
                for k=1:s1
                    J(sJ+1,n*(mm*ll+1)+1+nb_par+2*n+(k-1)*n+(1:n))=...
                        (step_cnd(j).v(:,k))';
                end;
                J(sJ+1,n*(mm*ll+1)+1+2*n+nb_par+s1*n+(1:s1))=step_cnd(j).lambda_v';
                
                s2=size(step_cnd(j).w,2);
                for k=1:s2
                    J(sJ+1,n*(mm*ll+1)+1+nb_par+2*n+s1*(n+1)+(k-1)*n+(1:n))=...
                        step_cnd(j).w(:,k)';
                end;
                J(sJ+1,n*(mm*ll+1)+1+2*n+nb_par+s1*(n+1)+s2*n+(1:s2))=...
                    step_cnd(j).lambda_w';
                J(sJ+1,n*(mm*ll+1)+1+2*n+nb_par+s1*(n+1)+s2*(n+1)+(1:s1))=...
                    (step_cnd(j).alpha)';
                res(sJ+1)=0;
            end;
            % add extra conditions
            if method.extra_condition
                [resi,condi]=sys_cond(point);
                for j=1:length(condi)
                    k=size(J,1);
                    res(k+1)=resi(j);
                    for o=1:ll*mm+1
                        J(k+1,(o-1)*n+1:o*n)=condi(j).profile(:,o)';
                    end;
                    J(k+1,n*mm*ll+n+1)=condi(j).period;
                    for o=1:length(free_par)
                        J(k+1,n*mm*ll+n+1+o)=condi(j).parameter(free_par(o));
                    end;
                    nb_par=length(free_par);
                    J(k+1,n*(mm*ll+1)+1+nb_par+(1:n))=condi(j).x1';
                    J(k+1,n*(mm*ll+1)+1+nb_par+n+(1:n))=condi(j).x2';
                    s1=length(condi(j).v(1,:));
                    for o=1:s1
                        J(k+1,n*(mm*ll+1)+1+nb_par+2*n+(o-1)*n+(1:n))=condi(j).v(:,o)';
                    end;
                    J(k+1,n*(mm*ll+1)+1+2*n+nb_par+s1*n+(1:s1))=condi(j).lambda_v';
                    s2=size(condi(j).w,2);
                    for o=1:s2
                        J(k+1,n*(mm*ll+1)+1+nb_par+2*n+s1*(n+1)+(o-1)*n+(1:n))=...
                            condi(j).w(:,o)';
                    end;
                    J(k+1,n*(mm*ll+1)+1+2*n+nb_par+s1*(n+1)+s2*n+(1:s2))=...
                        condi(j).lambda_w';
                    J(k+1,n*(mm*ll+1)+1+2*n+nb_par+s1*(n+1)+s2*(n+1)+(1:s1))=...
                        condi(j).alpha';
                end;
            end;
            % solve linear system
            if size(J,1)~=size(J,2)
                warning('p_correc:nonsquare','P_CORREC warning: use of nonsquare Jacobian.');
            end;
            dx=J\res;
            % apply non-parameter corrections
            for k=1:ll*mm+1
                point.profile(:,k)=point.profile(:,k)-real(dx((k-1)*n+1:k*n));
            end;
            point.period=point.period-real(dx(n*mm*ll+1+n));
            point.x1=point.x1-real(dx(n*mm*ll+1+n+(1:n)+length(free_par)));
            point.x2=point.x2-real(dx(n*mm*ll+1+2*n+(1:n)+length(free_par)));
            if ~isempty(point.v)
                s1=length(point.v(1,:));
            else
                s1=0;
            end;
            for k=1:s1
                point.v(:,k)=point.v(:,k) ...
                    -dx(n*mm*ll+n+1+length(free_par)+2*n+(k-1)*n+(1:n));
                point.lambda_v(k)=point.lambda_v(k) ...
                    -dx(n*mm*ll+n+1+length(free_par)+2*n+s1*n+k);
            end;
            if ~isempty(point.w)
                s2=length(point.w(1,:));
            else
                s2=0;
            end;
            for k=1:s2
                point.w(:,k)=point.w(:,k)-...
                    dx(n*mm*ll+n+1+length(free_par)+2*n+(k-1)*n+(1:n)+s1*(n+1));
                % complex conjugate here because the Jacobian is formulated in terms of the complex
                % conjugate of the residual... (NOT for the w's)
                point.lambda_w(k)=point.lambda_w(k)-conj(dx(n*mm*ll+n+1+length(free_par)+...
                    2*n+s1*(n+1)+s2*n+k));
            end;
            point.alpha=point.alpha-...
                dx(n*mm*ll+n+1+length(free_par)+2*n+length(point.v(1,:))*(n+1)+s2*(n+1)+(1:s1));
            %     % apply parameter corrections:
            for j=1:length(free_par)
                par(free_par(j))=par(free_par(j))-real(dx(p_start+j));
            end;
            
            % fill in parameters:
            point.parameter=par;
            

        otherwise
            p_from_x=str2func(['dde_',point.kind,'_from_x']);
            x_from_p=str2func(['dde_x_from_',point.kind]);
            kind_jac=str2func([point.kind,'_jac_combined']);
            jaccond=x_from_p([condi;step_cnd],free_par)';
            [J0,res0]=kind_jac(funcs,point,free_par,'reference',previous);
            res=[res0;resi;zeros(length(step_cnd),1)];
            J=[J0;jaccond];
            % solve linear system
            if size(J,1)~=size(J,2)
                warning('p_correc:nonsquare','P_CORREC warning: use of nonsquare Jacobian.');
            end;
            dx=J\res;
            x=x_from_p(point,free_par);
            point=p_from_x(x-dx,point,free_par);
            
    end
    
%     % apply parameter corrections:
%     for j=1:length(free_par)
%         par(free_par(j))=par(free_par(j))-real(dx(p_start+j));
%     end;
%     
%     % fill in parameters:
%     point.parameter=par;
%     
    % compute residual:
    norm_res=norm(res,'inf');
    if print_r
        fprintf('it=%d, res=%g\n',i, norm_res);
    end;
    
    % check for convergence
    if i==1
        r1=norm_res;
        r0=r1/2;
    else
        r0=r1;
        r1=norm_res;
    end;
    if r1<=conv_r || (i-1>nmon_iter && r1>=r0) || any(isnan(r1)) || any(isnan(dx))
        break;
    end;
    
end;

% compute final residual

n_res=norm(res,'inf');

success=(n_res<=method.minimal_accuracy) && all(~isnan(dx));

% recorrect with adapted mesh if necessary

if success && (strcmp(point.kind,'psol') || strcmp(point.kind,'hcli'))
    if size(point.mesh,2)
        ma=method.adapt_mesh_after_correct;
        if p_nr==1 || ( p_nr>1 && mod(p_nr,ma)==0 )
            % do not adapt mesh when p_nr=0
            % adapt mesh when p_nr=1
            % adapt mesh when p_nr>1 & mod(p_nr,ma)=0
            if success % adapt & correct again
                method2=method;
                method2.adapt_mesh_before_correct=1;
                method2.adapt_mesh_after_correct=0;
                [point,success]=p_correc(funcs,point,free_par,step_cnd,method2,2,previous);
            end
        end
    end
end

if isfield(method,'print_jacobian_condition') && exist('J','var')
    fprintf('norm(J^(-1))=%g\n',1/min(svd(J)));
end
end