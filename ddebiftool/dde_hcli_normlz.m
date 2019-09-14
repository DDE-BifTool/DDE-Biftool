function hcli=dde_hcli_normlz(hcli)
%% normalize hcli derived points
% (descended from |'hcli'| kind)
%
% $Id: dde_hcli_normlz.m 366 2019-07-14 21:52:54Z jansieber $
%%
for k=1:size(hcli.v,2)
    hcli.v(:,k)=hcli.v(:,k)/norm(hcli.v(:,k));
end
for k=1:size(hcli.w,2)
    hcli.w(:,k)=hcli.w(:,k)/norm(hcli.w(:,k));
end
if ~isempty(hcli.alpha)
    hcli.alpha=hcli.alpha/norm(hcli.alpha);
end
end
