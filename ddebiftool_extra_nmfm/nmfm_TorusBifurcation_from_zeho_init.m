function [psolbr,augmented] = nmfm_TorusBifurcation_from_zeho_init(funcs,zeho,radius,freepars,varargin)
%% Initialize branch for continuing the Neimark-Sacker curve
% emanating from the the zero-Hopf point (if at all).
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% (adapted by JS)
%
% $Id: nmfm_TorusBifurcation_from_zeho_init.m 369 2019-08-27 00:07:02Z jansieber $
%
if ~zeho.nmfm.transcritical
    %% generic case
    %% check condition for Neimark-Sacker bifurcation
    default={'fix_rotation',true};
    [options,pass_on]=dde_set_options(default,varargin,'pass_on');
    if real(zeho.nmfm.g011)*real(zeho.nmfm.g110)>0
        psolbr=repmat(df_brnch(funcs,freepars,'psol'),1,0);
        augmented=false;
        return
    end
    psoltemplate=p_topsol([],p_tohopf(funcs,zeho),...
        'radius',radius(1),pass_on{:});
    %% normal form coefficients & eigenvectors
    n=length(zeho.x);
    get_n=@(x)x(1:n);
    dev_eval=@(x)get_n(nmfm_dev_call(x,0));
    g111=real(zeho.nmfm.g111); % should be real anyway, but may be affected by roundoff
    g021=zeho.nmfm.g021;
    g200=zeho.nmfm.g200; % safely real
    g110=zeho.nmfm.g110;
    g011=real(zeho.nmfm.g011); % should be real anyway, but may be affected by roundoff
    
    K=zeho.nmfm.K;
    
    h000mu=[dev_eval(zeho.nmfm.h000mu(1)), dev_eval(zeho.nmfm.h000mu(2))]; % safely real
    h011  = real(dev_eval(zeho.nmfm.h011)); % should be real anyway, but may be affected by roundoff
    h020  = dev_eval(zeho.nmfm.h020);
    
    q0=zeho.nvec.q0;
    q1=zeho.nvec.q1;
    
    omega1=zeho.nmfm.omega1;
    omega2=zeho.nmfm.omega2;
    %% parameter expansion
    beta(1)=-g011;
    beta(2)=(2*(real(g110)-g200)*real(g021)+real(g110)*g111)/(2*g200);
    dpar=[zeho.parameter(freepars)', zeros(length(freepars),1), K*beta(:)];
    %% frequency expansion
    z0=-(2*real(g021)+g111)/(2*g200); % safely real
    d_om2=omega1*beta(1)+omega2*beta(2)+imag(g021)+imag(g110)*z0;
    d_om=[zeho.nvec.omega, 0, d_om2];
    %% profile expansion
    zn=zeros(size(zeho.x));
    profile={...
        [zeho.x, zn,   -g011*h000mu(:,1)+beta(2)*h000mu(:,2)+z0*q0+h011],...
        [zn,     2*q1],...
        [zn,     zn,   h020] };
    %% approximation to limit cycle
    psolbr=nmfm_psol_from_C2(funcs,profile,dpar,d_om,freepars,radius,psoltemplate);
    psolbr=replace_branch_pars(psolbr,psolbr.parameter.free,pass_on);
    augmented=false;
    if options.fix_rotation
        psolbr.tangent=@fix_rotation;
    end
else
    %% transcritical case: two Neimark-Sacker bifurcations curves
    % 
    
    %% check condition for Neimark-Sacker bifurcation
    default={'fix_rotation',true};
    [options,pass_on]=dde_set_options(default,varargin,'pass_on');
    if real(zeho.nmfm.g011)*real(zeho.nmfm.g200)<0
        psolbr=repmat(df_brnch(funcs,freepars,'psol'),1,0);
        augmented=false;
        return
    end
    psoltemplate=p_topsol([],p_tohopf(funcs,zeho),...
        'radius',radius(1),pass_on{:});
    %% normal form coefficients & eigenvectors
    n=length(zeho.x);
    get_n=@(x)x(1:n);
    dev_eval=@(x)get_n(nmfm_dev_call(x,0));
    g111=real(zeho.nmfm.g111); % should be real anyway, but may be affected by roundoff
    g021=zeho.nmfm.g021;
    g200=zeho.nmfm.g200; % safely real
    g110=zeho.nmfm.g110;
    g011=real(zeho.nmfm.g011); % should be real anyway, but may be affected by roundoff
    g210=zeho.nmfm.g210;
    g300=zeho.nmfm.g300;
    % parameter transformation
    K=zeho.nmfm.K;
    % centermanifold
%     h100mu=[dev_eval(zeho.nmfm.h100mu(1)), dev_eval(zeho.nmfm.h100mu(2))]; % safely real
%     h010mu=[dev_eval(zeho.nmfm.h010mu(1)), dev_eval(zeho.nmfm.h010mu(2))]; % safely real
    h011  = real(dev_eval(zeho.nmfm.h011)); % should be real anyway, but may be affected by roundoff
    h020  = dev_eval(zeho.nmfm.h020);
    h200  = dev_eval(zeho.nmfm.h200);
    h110  = dev_eval(zeho.nmfm.h200);
    % eigenvectors
    q0=zeho.nvec.q0;
    q1=zeho.nvec.q1;
    % coefficients of the period
    omega1=zeho.nmfm.omega1;
    omega2=zeho.nmfm.omega2;
    for i=[1,-1]
        z0(1)=i*sqrt(g011/g200);
        z0(2)=-(2*g200*real(g021)-g300*g011)/(2*g200^2);
        %% parameter expansion
        beta(1,1)=-i*2*sqrt(g011*g200);
        beta(1,2)=-(g300*g011/g200+g111);
        beta(2,1)=-real(g110)*z0(1);
        %beta(2,2)=-real(g110)*z0(1)^2-real(g210)*(z0(1)^2+z0(2))-real(g021);
        beta(2,2)=(-g011*real(g210)/g200-g300*g011*real(g110)/(2*g200^2)...
            + real(g110)*real(g021)/g200-real(g021));
        beta(1,2)=0; beta(2,2)=0;
        dpar=[zeho.parameter(freepars)', K*beta(:,1), K*beta(:,2)];
        %% frequency expansion
        d_om(1)=beta(1,1)*omega1+beta(2,1)*omega2+imag(g110)*z0(1);
        d_om(2)=beta(1,2)*omega1+beta(2,2)*omega2+imag(g110)*z0(2)+imag(g210)*z0(1)^2+imag(g021);
        d_om=[zeho.nvec.omega, d_om(1), d_om(2)];
        %% profile expansion
        zn=zeros(size(zeho.x));
%         profile={...
%             [zeho.x, z0(1)*q0, ...
%                 h011-2*g011*h100mu(:,1)+g011/(2*g200)*h200...
%                 -real(g021)/g200*q0-real(g110)*g011/g200*h100mu(:,2)],...
%             [zn,  2*q1, i*2*sqrt(g011*g200)*h010mu(:,1)...
%                         -2*z0(1)*h110+2*real(g110)*z0(1)*h010mu(:,2)],...
%             [zn,  zn, h020] };
        profile={...
            [zeho.x, z0(1)*q0], ... 
            [zn,  2*q1]};    
        %% approximation to limit cycle
        psolbr(1/2*i+3/2)=nmfm_psol_from_C2(funcs,profile,dpar,d_om,freepars,radius,psoltemplate);
        psolbr(1/2*i+3/2)=replace_branch_pars(psolbr(1/2*i+3/2),psolbr(1/2*i+3/2).parameter.free,pass_on);
    end
    if options.fix_rotation
      for i=[1,-1]
        psolbr(1/2*i+3/2).tangent=@fix_rotation;
      end
    end
    augmented=false;
end
end

function tangent=fix_rotation(p)
tangent=p_axpy(0,p(1),[]);
tangent.parameter(end-1)=1;
end
