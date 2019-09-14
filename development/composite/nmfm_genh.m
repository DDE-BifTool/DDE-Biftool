function newpoint = nmfm_genh(funcs, point, varargin)
%% Compute 2nd Lyapunov coefficient for generalized Hopf point
%
% $Id: nmfm_genh.m 186 2017-04-13 06:06:31Z mmbosschaert $
% 
%%
default={'nullpoint',[],'free_pars',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
kind = point.kind;
n = length(point.x);
newpoint = point;
if ~strcmp(kind,'genh')
    display(kind);
    error('NMFM_GENH: did not receive a generalized hopf point as argument.');
end
%% Select eigenvalue pair
omega = point.omega;
if isempty(omega) || omega == 0
    fprintf('NMFM_GENH: omega is empty or zero, returning L1 = NaN.\n');
    newpoint.nmfm.L2 = NaN;
    return;
end
lambda0 = 1i*omega;
%% Construct characteristic matrix
Delta=@(z,order)ch_matrix(funcs,point.x,point.parameter,z,'deri',order);
Delta2 = Delta(2*lambda0,0);
Delta3 = Delta(3*lambda0,0);
Delta0 = Delta(0,0);
DDelta0 = Delta(0,1);
DDelta2 = Delta(2*lambda0,1);
D2Delta1 = Delta(lambda0,2);
%% Compute nullvectors
[p,q]=nmfm_nullvector(funcs,point,lambda0,'nullpoint',options.nullpoint,pass_on{:});
%% abbreviate zero parameter vector and define/abbreviate derivatives
par0=point.parameter(:)*0;
F=nmfm_deriv_define(funcs,point,'free_pars',options.free_pars,pass_on{:});
%% Implement normal form computations
phi = nmfm_dev_fun([q;par0],'lambda',lambda0);
phibar = nmfm_dev_conj(phi);

D2_by_B_phi_phi=Delta2\F.B(phi,phi);
D0_by_B_phi_phibar=Delta0\F.B(phi,phibar);
h20 = nmfm_dev_fun([D2_by_B_phi_phi;par0],'lambda',2*lambda0);
h20bar=nmfm_dev_conj(h20);
h11 = nmfm_dev_fun([D0_by_B_phi_phibar;par0]);

c1 = (1/2)*p*(F.B(phibar,h20) + 2*F.B(phi,h11) + F.C(phi,phi,phibar));

h30=nmfm_dev_fun([Delta3\( 3*F.B(phi,h20) + F.C(phi,phi,phi) );par0],'lambda',3*lambda0);
h21=nmfm_dev_fun( [(-2*c1)*[((1/2)*p*D2Delta1*q)*q,-q];par0,par0], 'lambda',lambda0*[1,1],'t',[0,1]);
h21bar=nmfm_dev_conj(h21);

h31_v1=Delta2\(...
    F.B(phibar,h30)+3*F.B(phi,h21)+3*F.B(h20,h11)+...
    3*F.C(phi,phibar,h20)+3*F.C(phi,phi,h11)+...
    F.D(phi,phi,phi,phibar));
h31_v20=6*c1*(Delta2\(DDelta2-eye(n)))*D2_by_B_phi_phi;
h31_v21=-6*c1*D2_by_B_phi_phi;
h31=nmfm_dev_fun([h31_v1,h31_v20,h31_v21;par0(:,[1,1,1])],'lambda',2*lambda0*[1,1,1],'t',[0,0,1]);

h22 = nmfm_dev_fun( [Delta0\(...
    2*F.B(phibar,h21) + 2*F.B(h11,h11) + 2*F.B(phi,h21bar) + F.B(h20,h20bar) + ...
    F.C(phibar,phibar,h20) + F.C(phi,phi,h20bar) + 4*F.C(phi,phibar,h11) +...
    F.D(phi,phi,phibar,phibar)); par0] );

c2 = (1/12)*p*(6*F.B(h11,h21) + 3*F.B(h21bar,h20) + F.B(h20bar,h30) + 3*F.B(phi,h22) + 2*F.B(phibar,h31) + ...
    6*F.C(phibar,h20,h11) + 6*F.C(phi,h11,h11) + 3*F.C(phi,h20,h20bar) + 6*F.C(phi,phibar,h21) + ...
    3*F.C(phi,phi,h21bar) + F.C(phibar,phibar,h30) + ...
    6*F.D(phi,phi,phibar,h11) + 3*F.D(phi,phibar,phibar,h20) + ...
    F.D(phi,phi,phi,h20bar) + ...
    F.E(phi,phi,phi,phibar,phibar));

L2 = real(c2)/(omega);
%% Attach result to point structure
if ~isfield(newpoint.nmfm,'L1') || isempty(newpoint.nmfm.L1) || isnan(newpoint.nmfm.L1)
    newpoint.nmfm.L1 = real(c1)/omega;
end

newpoint.nmfm.L2 = L2;
newpoint.nvec.p = p;
newpoint.nvec.q = q;

%% parameter-related normal form coefficients
if ~isempty(options.free_pars)
    vp=eye(2);
    h00mu=[nmfm_dev_fun([Delta0\F.J1*vp(:,1);par0]),...
        nmfm_dev_fun([Delta0\F.J1*vp(:,2);par0])];
    
    gamma1mu=arrayfun(@(i)p*(F.A1(phi,vp(:,i))+F.B(phi,h00mu(i))),1:2);
    
    h10mu=arrayfun(@(i)nmfm_binv(funcs,point,lambda0,q,p,...
        F.A1(phi,vp(:,i))+F.B(phi,h00mu(i)),-gamma1mu(i)),1:2);
    h01mu=[nmfm_dev_conj(h10mu(1)),nmfm_dev_conj(h10mu(2))];
    
    I=eye(n);
    h20_0=nmfm_dev_call(h20,0);
    h20_0=h20_0(1:n);
    h11_0=nmfm_dev_call(h11,0);
    h11_0=h11_0(1:n);
    for i=2:-1:1
        h20mu_exp(:,i) = DDelta2\(...
            -2*gamma1mu(i)*h20_0 + F.A1(h20,vp(:,i)) + 2*F.B(phi,h10mu(i)) +...
            F.B(h20,h00mu(i)) + F.B1(phi,phi,vp(:,i)) + F.C(phi,phi,h00mu(i)));
        h20mu_0(:,i) = -2*gamma1mu(i)*(DDelta2-I)*(Delta2\F.B(phi,phi));
        h20mu_1(:,i) = -2*gamma1mu(i)*DDelta2\F.B(phi,phi);
        h20mu(i) = nmfm_dev_fun([h20mu_exp(:,i),h20mu_0(:,i),h20mu_1(:,i);par0(:,[1,1,1])],...
            'lambda',[2*lambda0,0,0],'t',[0,0,1]);
        %%
        h11mu_0(:,i) = Delta0\(...
            -2*real(gamma1mu(i))*h11_0 + F.B(phibar,h10mu(i)) + F.B1(phi,phibar,vp(:,i))+...
            F.C(phi,phibar,h00mu(i)) + F.A1(h11,vp(:,i)) + F.B(phi,h01mu(i)) + F.B(h11,h00mu(i)))-...
            2*real(gamma1mu(i))*(DDelta0-I)*(Delta0\F.B(phi,phi));
        h11mu_1(:,i) = 2*real(gamma1mu(i))*Delta0\F.B(phibar,phi);
        h11mu(i) = nmfm_dev_fun([h11mu_0(:,i),h11mu_1(:,i);par0(:,[1,1])],...
            'lambda',[0,0],'t',[0,1]);
        %%
        gamma2mu(i) = ...
            1/2*p*(F.A1(h21,vp(:,i))+F.B(phibar,h20mu(i))+2*F.B(phi,h11mu(i))...
            +F.B(h21,h00mu(i))+F.B(h20,h01mu(i))+2*F.B(h11,h10mu(i))...
            +F.B1(h20,phibar,vp(:,i))+2*F.B1(phi,h11,vp(:,i))+2*F.C(phi,phibar,h10mu(i))...
            +F.C(h20,phibar,h00mu(i))+F.C(phi,phi,h01mu(i))+2*F.C(phi,h11,h00mu(i))...
            +F.C1(phi,phi,phibar,vp(:,i))+F.D(phi,phi,phibar,h00mu(i)));
    end
    % b1
    K01=real([gamma1mu(1), gamma1mu(2); gamma2mu(1), gamma2mu(2)])\[0;1];
    b12=imag([gamma1mu(1),gamma1mu(2)]*K01);
    
    newpoint.nmfm.b12=b12;
    newpoint.nmfm.c1=c1;
    newpoint.nmfm.c2=c2;
    
    newpoint.nmfm.PHI=phi;
    newpoint.nmfm.K01=K01;
    newpoint.nmfm.H1001=h10mu(1);
    newpoint.nmfm.H0001=h00mu(2);
    newpoint.nmfm.H11=h11;
    newpoint.nmfm.H20=h20;
    newpoint.nmfm.H30=h30;
    newpoint.nmfm.H21=h21;
end
end
