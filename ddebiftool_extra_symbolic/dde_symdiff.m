function [df,xxdev,pardev]=dde_symdiff(fs,xxs,ps,varargin)
%% directional derviatives of fs wrt xxs and ps
%
%% Inputs:
%
% * fs: m x 1 array of symbolic expressions using xxs as states and ps as
% parameters
% * xxs: n x ntau+1 array of current and delayed states
% * ps: 1 x np array of parameter symbols
%
%% Optional inputs (name-value pairs) 
%
% * maxorder: 5, 
% * append='_dev': names used for deviation directions, make sure these do
% not clash with original names in xxs and ps.
% * deviation_name: default |h_devsmall| deviation variable (will be set to
% 0 after differentiation)
%
%% Outputs:
%
% * df: cell array of length maxorder where |df{k}| contains the the kth
% directional derivative.
% $Id: dde_symdiff.m 309 2018-10-28 19:02:42Z jansieber $
%%
default={'maxorder',5,'dev_append','_d','deviation_name','h_devsmall','directional_derivative',true};
options=dde_set_options(default,varargin,'pass_on');
[n,ntausp1]=size(xxs);
npar=length(ps);
%% extract variable names
xnames=dde_names_from_sym(xxs);
pnames=dde_names_from_sym(ps);
xpnames=[xnames(:).',pnames(:).'];
%% define deviations
if options.directional_derivative
    rep=1;
else
    rep=options.maxorder;
end
hname=options.deviation_name;
xdevnames=dde_name_append(xnames,options.dev_append,rep);
pdevnames=dde_name_append(pnames,options.dev_append,rep);
devnames=[xdevnames(:).',pdevnames(:).'];
if ~isempty(intersect(xpnames,devnames)) || ismember(hname,[xpnames,devnames])
    error('symdiff:names',...
        'symdiff: name clash: choose different option for append or deviation_name\n%s\n%s\n%s\n',...
        evalc('celldisp(xpnames,''variables'')'),...
        evalc('celldisp(xpnames,''deviations'')'),...
        evalc('celldisp(xpnames,''hname'')'));
end
for i=rep:-1:1
    xxdev{i}=reshape(cell2sym(xdevnames(:,:,i)),[n,ntausp1]);
    pardev{i}=reshape(cell2sym(pdevnames(:,:,i)),[1,npar]);
end
hdev=sym(options.deviation_name);
deriv=fs;
df=cell(1,options.maxorder);
%% differentiate i times
for i=1:options.maxorder
    if options.directional_derivative
        irep=1;
    else
        irep=i;
    end
    %% add deviations to xx and par
    fdev0=subs(deriv,xxs,xxs+hdev*xxdev{irep});
    fdev=subs(fdev0,ps,ps+hdev*pardev{irep});
    deriv=diff(fdev,hdev);
    df{i}=subs(deriv,hdev,0);
    deriv=df{i};
end
end
%% append letter to cell array of names
function xnap=dde_name_append(xn,app,rep)
xnap=cell([numel(xn),rep]);
for k=1:rep
    rapp=repmat(app,1,k);
    for i=1:numel(xn)
        xnap{i,k}=[xn{i},rapp];
    end
end
xnap=reshape(xnap,[size(xn),rep]);
end