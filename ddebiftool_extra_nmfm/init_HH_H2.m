function hopf_br = init_HH_H2(funcs,hoho,eps)
%% Initialize branch for continuing the second Hopf curve
% emanating from the Hopf-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_HH_H2.m 309 2018-10-28 19:02:42Z jansieber $
%%

freepars=get_free_pars();
K01=hoho.nmfm.K01;

hopf1=p_tohopf(funcs,hoho,hoho.omega2);
hopf2=hopf1;

hopf2.parameter(freepars)=hopf2.parameter(freepars)+eps*K01';

method=df_mthod(funcs,'hopf');
[hopf2,s]=p_correc(funcs,hopf2,freepars(2),[],method.point);
if ~s
    warning('Could not correct second fold point.')
end

hopf_br=df_brnch(funcs,freepars,'hopf');
hopf_br.point=hopf1;
hopf_br.point(2)=hopf2;

end
