function hopf_br = init_HH_H1(funcs,hoho,eps)
%% Initialize branch for continuing the first Hopf curve
% emanating from the Hopf-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_HH_H1.m 309 2018-10-28 19:02:42Z jansieber $
%%

freepars=get_free_pars();
K10=hoho.nmfm.K10;

hopf1=p_tohopf(funcs,hoho,hoho.omega1);
hopf2=hopf1;

hopf2.parameter(freepars)=hopf2.parameter(freepars)+eps*K10';

method=df_mthod(funcs,'hopf');
[hopf2,s]=p_correc(funcs,hopf2,freepars(2),[],method.point);
if ~s
    warning('Could not correct second Hopf point.')
end

hopf_br=df_brnch(funcs,freepars,'hopf');
hopf_br.point=hopf1;
hopf_br.point(2)=hopf2;

end
