function [pdevs,factors]=nmfm_pol_from_devs(devs,order,groups)
%% decompose mixed derivative into linear combination of directional derivatives
%
% using the polarization identity
%
%% Inputs
%
% * devs: 1 x npol array of history functions (struct's)
% * order: requested order of derivative
% * groups (cell array): decomposition of 1:npol into groups such that devs
% in each group are equal (e.g., D3f(x)[p,p,q] would have groups={[1,2],3}
% while D3f(x)[p,q,r] should have groups={1,2,3})
%
%% Outputs
%
% * pdevs: deviations in which one has to compute the directional
% derivatives
% * factors: coefficients of linear combination of pdevs to obtain mixed
% derivative with devs.
% 
% $Id$
%#ok<*AGROW>
%% for test purposes (devs can be real or complex numbers or dev functions)
if isnumeric(devs)
    ax=@(c,dev)sum(c.*dev);
    toconj=@(dev)conj(dev);
    isreal=@(dev)imag(dev)==0;
else
    ax=@nmfm_dev_ax;
    toconj=@nmfm_dev_conj;
    isreal=@(dev)all(imag(dev.lambda)==0);
end
%% determine polarization sums and factors
[cmats,cfactors]=dde_polarization_coeffs(order,groups);
%% pick one dev from each group of equals
devs=devs(cellfun(@(x)x(1),groups));
%% collect deviations to be applied to point (still possibly complex)
npol=length(cfactors);
for i=npol:-1:1
    cdevs(i)=ax(cmats(:,i).',devs);
end
%% for all complex deviations, split up into real deviations
[mat_cexp,fac_cexp]=nmfm_pol_real_from_complex(order);
pdevs=[];
factors=[];
for i=1:npol
    if isreal(cdevs(i))
        pdevs=[pdevs,cdevs(i)];
        factors=[factors,cfactors(i)];
    else %% complex cdevs(i)
        rpart=ax(0.5*[1,1],[cdevs(i),toconj(cdevs(i))]);
        ipart=ax(0.5*[-1i,1i],[cdevs(i),toconj(cdevs(i))]);
        for k=length(fac_cexp):-1:1
            rdevs(k)=ax(mat_cexp(:,k).',[rpart,ipart]);
        end
        pdevs=[pdevs,rdevs];
        factors=[factors,cfactors(i)*fac_cexp];
    end
end
end
