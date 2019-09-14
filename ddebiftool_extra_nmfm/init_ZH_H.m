function hopf_br = init_ZH_H(funcs,zeho,eps)
%% Initialize branch for continuing the Hopf curve
% emanating from the fold-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_ZH_H.m 314 2019-01-24 14:28:23Z mmbosschaert $
%%
freepars=get_free_pars();
K01=zeho.nmfm.K(:,2);

hopf1=p_tohopf(funcs,zeho);
hopf2=hopf1;

hopf2.parameter(freepars)=hopf2.parameter(freepars)+eps*K01';

method=df_mthod(funcs,'hopf');
[hopf2,s]=p_correc(funcs,hopf2,freepars(2),[],method.point);
if ~s
    warning('Could not correct second Hopf point.')
end

hopf_br=df_brnch(funcs,freepars,'hopf');
hopf_br.point=hopf1;
hopf_br.point(2)=hopf2;

end
