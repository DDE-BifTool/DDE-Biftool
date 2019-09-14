function point=dde_TorusBifurcation_postprocess(point,data)
%% renormalize eigenvector to be as close as possible to reference
%
% $Id: dde_TorusBifurcation_postprocess.m 369 2019-08-27 00:07:02Z jansieber $
%%
prev=data.previous;
n=size(prev.profile,1)/3;
ivre=n+(1:n);
ivim=2*n+(1:n);
prev_re=setfield(prev,'profile',prev.profile(ivre,:)); %#ok<*SFLD>
prev_im=setfield(prev,'profile',prev.profile(ivim,:));
vre=setfield(point,'profile',point.profile(ivre,:));
vim=setfield(point,'profile',point.profile(ivim,:));
sc=@(p1,p2)dde_coll_profile_dot(p1,p2);
nrm=sqrt(sc(vre,vre)+sc(vim,vim));
vre.profile=vre.profile/nrm;
vim.profile=vim.profile/nrm;
c=sc(vre,prev_re)+sc(vim,prev_im)+1i*(sc(vre,prev_im)-sc(vim,prev_re));
phi=angle(c);
vc=(vre.profile+1i*vim.profile)*exp(1i*phi);
point.profile(ivre,:)=real(vc);
point.profile(ivim,:)=imag(vc);
end

