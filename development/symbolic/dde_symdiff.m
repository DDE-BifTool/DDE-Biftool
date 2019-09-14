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
% $Id: dde_symdiff.m 170 2017-03-05 03:39:50Z jansieber $
%%
default={'maxorder',5,'dev_append','_dev','deviation_name','h_devsmall'};
options=dde_set_options(default,varargin,'pass_on');
[n,ntausp1]=size(xxs);
npar=length(ps);
fs=fs(:);
%% define deviations
xnames=arrayfun(@(x)char(x),xxs(:).','uniformoutput',false);
pnames=arrayfun(@(x)char(x),ps,'uniformoutput',false);
xdevnames=cellfun(@(x)[x,options.dev_append],xnames,'uniformoutput',false);
pdevnames=cellfun(@(x)[x,options.dev_append],pnames,'uniformoutput',false);
xpnames=[xnames(:).',pnames(:).'];
devnames=[xdevnames(:).',pdevnames(:).'];
hname=options.deviation_name;
if ~isempty(intersect(xpnames,devnames)) || ismember(hname,[xpnames,devnames])
    error('symdiff:names',...
        'symdiff: name clash: choose different option for append or deviation_name\n%s\n%s\n%s\n',...
        evalc('celldisp(xpnames,''variables'')'),...
        evalc('celldisp(xpnames,''deviations'')'),...
        evalc('celldisp(xpnames,''hname'')'));
end
xxdev=reshape(sym(xdevnames),[n,ntausp1]);
pardev=reshape(sym(pdevnames),[1,npar]);
hdev=sym(options.deviation_name);
%% add deviations to xx and par
fdev0=subs(fs,xxs,xxs+hdev*xxdev);
fdev=subs(fdev0,ps,ps+hdev*pardev);
%% differentiate i times
df=cell(1,options.maxorder);
for i=1:options.maxorder
    hi=repmat({hdev},1,i);
    deriv=diff(fdev,hi{:});
    df{i}=subs(deriv,hdev,0);
end
end