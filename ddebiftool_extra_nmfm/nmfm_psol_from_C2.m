function psolbr=nmfm_psol_from_C2(funcs,profile,param,omega,freepars,radius,psoltemplate)
%% create periodic orbit from expansion coefficients
%
% profile is cell array of length nf+1 if complex amplitudes of
% exp(1i*k*t/(2pi)) for k from 0 up to nf are needed. Each element of
% profile is a n x (order+1) array x where x(:,j) is the expansion
% coefficient of this harmonic amplitude in front of radius^(j-1)
%
% $Id: nmfm_psol_from_C2.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
maxfreqp1=length(profile);
psi=2*pi*psoltemplate.mesh;
r=ones(maxfreqp1,length(psi));
for k=maxfreqp1-1:-1:1
    r(k+1,:)=exp(1i*psi*k);
end
psol=repmat(psoltemplate,1,length(radius));
for nsol=1:length(radius)
    xc=zeros(size(profile{1},1),length(psoltemplate.mesh));
    par=zeros(size(param,1),1);
    om=0;
    for k=1:maxfreqp1
        for j=1:size(profile{k},2)
            xc = xc + profile{k}(:,j)*r(k,:)*radius(nsol)^(j-1);
        end
    end
    for j=1:size(param,2)
        par=par+param(:,j)*radius(nsol)^(j-1);
    end
    for j=1:length(omega)
        om=om+omega(j)*radius(nsol)^(j-1);
    end
    psol(nsol).parameter(freepars)=par';
    psol(nsol).profile=real(xc);
    psol(nsol).period=2*pi/om;    
end
psolbr=df_brnch(funcs,freepars,'psol');
psolbr.point=psol;
end
