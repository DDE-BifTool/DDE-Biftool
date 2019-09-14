function hopf_br = init_HT_H1(funcs,zeho,eps)
%% Initialize branch for continuing the first Hopf curve
% emanating from the Hopf-Transcritical point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_HT_H1.m 309 2018-10-28 19:02:42Z jansieber $
%%

freepars=get_free_pars();
K10=zeho.nmfm.K10;
K01=zeho.nmfm.K01;

hopf1=p_tohopf(funcs,zeho);
hopf2=hopf1;

beta1=eps;
beta2=0;
pm=K10*beta1+K01*beta2;
hopf2.parameter(freepars)=hopf2.parameter(freepars)+pm';

method=df_mthod(funcs,'hopf');
[hopf2,s]=p_correc(funcs,hopf2,freepars(2),[],method.point);
if ~s
    warning('Could not correct second Hopf point.')
end

hopf_br=df_brnch(funcs,freepars,'hopf');
hopf_br.point=hopf1;
hopf_br.point(2)=hopf2;

end
