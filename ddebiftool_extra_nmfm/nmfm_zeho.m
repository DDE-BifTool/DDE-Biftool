function zeho = nmfm_zeho(funcs, point, varargin)
%% Compute normal form for zero-Hopf interaction
%
%  $Id: nmfm_zeho.m 317 2019-01-30 17:12:33Z mmbosschaert $
%
%%
default={'nullpoint',[],'print',0,'free_pars',[],'transcritical',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
zeho = point;
if ~strcmp(zeho.kind,'zeho')
    error('NMFM_ZEHO: did not receive a zero hopf point, but %s, as argument.',zeho.kind);
end

%% Select eigenvalue pair
if isempty(zeho.nvec.omega) || zeho.nvec.omega == 0
    warning('NMFM_ZEHO:omega',...
        'NMFM_ZEHO: omega is empty or zero, unable to compute normal form.');
    return
end

lambda0 = 0; % We could try to take the actual eigenvalue approximately equal to 0
lambda1 = 1i*zeho.nvec.omega;

%% Compute nullvectors
par=point.parameter;
Delta=@(lambda)ch_matrix(funcs,point.x,par,lambda);
dDelta=@(lambda)ch_matrix(funcs,point.x,par,lambda,'deri',1);
binv=@(lambda,q,p,zeta,kappa)nmfm_binv(funcs,point,lambda, q, p, zeta, kappa);
[p0,q0]=nmfm_nullvector(funcs,point,lambda0,pass_on{:});
[p1,q1]=nmfm_nullvector(funcs,point,lambda1,pass_on{:});
%% abbreviate derivatives, characteristic matrix and bordered inverse
F=nmfm_deriv_define(funcs,point,...
    'free_pars',options.free_pars,'print',options.print,pass_on{:});
par0=par(:)*0;
dev0=@(v)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))]);
devl=@(v,lambda)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',lambda);
devlt=@(v,lambda,t)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',lambda,'t',t);

%% Eigenfunctions
phi0 = dev0(q0);
phi1 = devl(q1,lambda1);
phi1bar = nmfm_dev_conj(phi1);

%% Some often used derivatives
D2F00 = F.B(phi0,phi0);
D2F01 = F.B(phi0,phi1);
D2F11b =real(F.B(phi1,phi1bar)); % should be real but unsafe

%% Normal form coefficients
%% Quadratic
g200 = (1/2)*p0*D2F00;
g110 = p1*D2F01;
g011 = p0*D2F11b;

%% Center manifold
h200 = binv(lambda0,q0,p0,D2F00,-p0*D2F00);
h020 = devl( Delta(2*lambda1)\F.B(phi1,phi1),2*lambda1);
h110 = binv(lambda1,q1,p1,D2F01,-p1*D2F01);
h110bar=nmfm_dev_conj(h110);
h011 = binv(lambda0,q0,p0,D2F11b,-p0*D2F11b);

%% Cubic
g300 = (1/6)*p0*( 3*F.B(phi0,h200) + F.C(phi0,phi0,phi0) );
g111 = p0*( F.B(phi0,h011) + F.B(phi1bar,h110) + F.B(phi1,h110bar) + F.C(phi0,phi1,phi1bar));
g210 = (1/2)*p1*( F.B(phi1,h200)    + 2*F.B(phi0,h110) + F.C(phi0,phi0,phi1) );
g021 = (1/2)*p1*( F.B(phi1bar,h020) + 2*F.B(phi1,h011) + F.C(phi1,phi1,phi1bar) );

%% Gavrilov
b = g200;
c = g011;
d = g110 - lambda1*g300/g200;
e = real(g210 + g110*(real(g021)/g011 - 3*g300/(2*g200) + g111/(2*g011)) - g021*g200/g011);

s = b*c;
theta0 = real(g110)/g200;
if options.print>0
    fprintf('s = %.10f, theta = %.10f.\n', s, theta0);
end
%% Pass on normal form coefficients
zeho.nmfm.g200 = g200;
zeho.nmfm.g110 = g110;
zeho.nmfm.g011 = g011;
zeho.nmfm.g300 = g300;
zeho.nmfm.g111 = real(g111);
zeho.nmfm.g210 = g210;
zeho.nmfm.g021 = g021;
zeho.nmfm.b = b;
zeho.nmfm.c = c;
zeho.nmfm.d = d;
zeho.nmfm.e = e;
zeho.nmfm.s = s;
zeho.nmfm.theta = theta0;

zeho.nvec.p0 = p0;
zeho.nvec.q0 = q0;
zeho.nvec.p1 = p1;
zeho.nvec.q1 = q1;

%% parameter-related coefficients
if isempty(options.free_pars)
    return
end

zeho.nmfm.transcritical=options.transcritical;
zeho.nmfm.h200=h200;
zeho.nmfm.h011=h011;
zeho.nmfm.h020=h020;
zeho.nmfm.h110=h110;

if ~zeho.nmfm.transcritical % generic case
        
    D = [Delta(lambda0) q0 ; p0, 0];
    get_n=@(x)x(1:length(q0));
    
    gg = (p0*F.J1)';
    s1 = gg/(gg'*gg); 
    s2 = [gg(2) ; -gg(1)];

    r1=dev0(get_n(D\[F.J1*s1 ;0]));
    r2=dev0(get_n(D\[F.J1*s2 ;0]));
    r3=devlt([get_n(D\[dDelta(lambda0)*q0 ;0]),-q0],[0,0],[0,1]);    
    
    LL=[p0*(F.A1(phi0,s2) + F.B(phi0,r2)), p0*F.B(phi0,phi0);
        p1*(F.A1(phi1,s2) + F.B(phi1,r2)), p1*F.B(phi1,phi0)];
        
    RR = -[ p0*(F.A1(phi0,s1) + F.B(phi0,r1) - F.B(phi0,r3));
            p1*(F.A1(phi1,s1) + F.B(phi1,r1) - F.B(phi1,r3))];
        
    dd = real(LL)\[real(RR),[0;1]];
    K  = ([[1;0], dd(1,:)']*[s1,s2]')';
    
    % K10 = s1 + dd(1,1)*s2
    % K01 = dd(1,2)*s2
    
    h000mu(1)=nmfm_dev_ax( [1;dd(:,1);-1], [r1,r2,phi0,r3] );
    h000mu(2)=nmfm_dev_ax(  dd(:,2),       [r2,phi0] );    
    % coefficients omega1 and omega2
    omega1=imag(p1*F.B(phi1,h000mu(1)) + p1*F.A1(phi1,K(:,1)));
    omega2=imag(p1*F.B(phi1,h000mu(2)) + p1*F.A1(phi1,K(:,2)));
        
    zeho.nmfm.K=K;
    zeho.nmfm.h000mu=h000mu;
    zeho.nmfm.omega1=omega1;
    zeho.nmfm.omega2=omega2;
else % transcritical case
    vp=eye(2);
    LL=[p0*F.A1(phi0,vp(:,1)),p0*F.A1(phi0,vp(:,2));
        p1*F.A1(phi1,vp(:,1)),p1*F.A1(phi1,vp(:,2))];
    K=real(LL)\vp;
    
    % coefficients omega1 and omega2
    omega1=imag(p1*F.A1(phi1,K(:,1)));
    omega2=imag(p1*F.A1(phi1,K(:,2)));
    
%     h100mu(1)=binv(lambda0,q0,p0,F.A1(phi0,K(:,1)),-1);
%     h100mu(2)=binv(lambda0,q0,p0,F.A1(phi0,K(:,2)), 0);
%     h010mu(1)=binv(lambda0,q0,p0,F.A1(phi1,K(:,1)), -1i*omega1);
%     h010mu(2)=binv(lambda0,q0,p0,F.A1(phi1,K(:,2)),-(1+1i*omega2));
    
%     h002mu=devl(Delta(-2*lambda1)\F.B(phi1bar,phi1bar), -2*lambda1);
    %% return coefficients
    zeho.nmfm.K=K;
    
%     zeho.nmfm.h100mu=h100mu;
%     zeho.nmfm.h010mu=h010mu;
%     zeho.nmfm.h002mu=h002mu;
    
    zeho.nmfm.omega1=omega1;
    zeho.nmfm.omega2=omega2;
end
end
