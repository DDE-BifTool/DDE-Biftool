function [hopf_br,augmented] = nmfm_hopf_from_hoho_init(funcs,hoho,radius,freepars,varargin)
%% Initialize branch for continuing Hopf curve
% emanating from the Hopf-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
%
% (adapted by JS)
%
% $Id: nmfm_hopf_from_hoho_init.m 315 2019-01-29 19:42:21Z jansieber $
%%
% default={'freqs',hoho.omega};
% options=dde_set_options(default,varargin,'pass_on');
% [hopftemplate,ind]=hoho_tohopf(funcs,hoho,options.freqs);
% ibr=[ind,3-ind];
% for k=length(ibr):-1:1
%     Kvec=hoho.nmfm.K(:,ibr(k));
%     Kvec=Kvec/norm(Kvec);
%     hopf=repmat(hopftemplate,1,length(radius));
%     for i=1:length(radius)
%         hopf(i).parameter(freepars)=hoho.parameter(freepars)+radius(i)*Kvec';
%     end
%     hopf_br(k)=df_brnch(funcs,freepars,'hopf');
%     hopf_br(k).point=hopf;
% end
% augmented=true;
for k=2:-1:1
    default={'freqs',hoho.nvec.omega(k)};
    options=dde_set_options(default,varargin,'pass_on');
    [hopftemplate,ind]=hoho_tohopf(funcs,hoho,options.freqs);
    Kvec=hoho.nmfm.K(:,k);
    Kvec=Kvec/norm(Kvec);
    hopf=repmat(hopftemplate,1,length(radius));
    for i=1:length(radius)
        hopf(i).parameter(freepars)=hoho.parameter(freepars)+radius(i)*Kvec';
    end
    hopf_br(k)=df_brnch(funcs,freepars,'hopf');
    hopf_br(k).point=hopf;
    omega(k)=hopf_br(k).point(end).omega;
end
[dum,ix]=sort(abs(omega-hoho.omega),'descend');
hopf_br=hopf_br(ix);
augmented=true;
end
