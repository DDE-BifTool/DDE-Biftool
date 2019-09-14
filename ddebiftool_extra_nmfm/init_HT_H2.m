function hopf_br = init_HT_H2(funcs,ht,eps)
%% Initialize branch for continuing the second Hopf curve
% emanating from the Hopf-Transcritical point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_HT_H2.m 309 2018-10-28 19:02:42Z jansieber $
%%

freepars=get_free_pars();
K10=ht.nmfm.K10;
K01=ht.nmfm.K01;
g110=ht.nmfm.g110;
g200=ht.nmfm.g200;

hopf1=p_tohopf(funcs,ht);
hopf2=hopf1;

beta1=eps;
beta2=real(g110)/g200*eps;
pm=K10*beta1+K01*beta2;
hopf2.parameter(freepars)=hopf2.parameter(freepars)+pm';
hopf2.x=hopf2.x-beta2/real(g110)*ht.nvec.q0;

method=df_mthod(funcs,'hopf');
[hopf2,s]=p_correc(funcs,hopf2,freepars(2),[],method.point);
if ~s
    warning('Could not correct second Hopf point.')
end

hopf_br=df_brnch(funcs,freepars,'hopf');
hopf_br.point=hopf1;
hopf_br.point(2)=hopf2;

end
