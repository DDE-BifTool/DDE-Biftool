function point=dde_hopf_postprocess(point,data)
%% change complex phase of Hopf eigenvector after newton iteration
% making it as close to the reference as possible
%
% $Id: dde_hopf_postprocess.m 308 2018-10-28 15:08:12Z jansieber $
%%
vref=data.previous.v;
v=point.v;
phi=angle(v'*vref);
point.v=v*exp(1i*phi);
end
