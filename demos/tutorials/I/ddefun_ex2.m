function  f = ddefun_ex2(xx,par)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: ddefun_ex2.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
    f = par(1)*xx(1,2)/(1+xx(1,2)^par(3))-par(2)*xx(1,1);
    % par = [beta, gamma, n];
end
