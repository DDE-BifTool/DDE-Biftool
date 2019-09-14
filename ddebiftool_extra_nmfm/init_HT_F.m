function fold_br = init_HT_F(funcs,zeho,eps)
%% Initialize branch for continuing the transcritical curve
% emanating from the Hopf-Transcritical point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_HT_F.m 314 2019-01-24 14:28:23Z mmbosschaert $
%%

freepars=get_free_pars();
K10=zeho.nmfm.K(:,1);
K01=zeho.nmfm.K(:,2);

fold1=p_tofold(funcs,zeho);
fold2=fold1;

% this needs to be fixed, the predictor should be
% fold2.parameter(freepars)=fold2.parameter(freepars)+eps*K01';
fold2.parameter(freepars)=fold2.parameter(freepars)+eps*K10';
%%
method=df_mthod(funcs,'fold');
[fold3,s]=p_correc(funcs,fold2,freepars(1),[],method.point);
if ~s
    warning('Could not correct second fold point.')
    fold3=fold2;
end

fold_br=df_brnch(funcs,freepars,'fold');
fold_br.point=fold1;
fold_br.point(2)=fold3;

end
