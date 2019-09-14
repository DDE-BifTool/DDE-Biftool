function [dfeq,xdevnames]=dde_sdmf_symdiff(fs,taus,xxs,varargin)
%% directional derviatives of f(xt) as defined by fs and taus wrt xxs
%
%% Inputs:
%
% * fs: m x 1 array of symbolic expressions using xxs as states
% * taus: ntau x 1 array of symbolic expressions using xxs as states
% * xxs: n x ntau+1 array of current and delayed states
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
% * dfeqs: cell array of length maxorder where |df{k}| contains the the kth
% directional derivative.
% * xdevf: symbolic expressions for deviations (cell array of symfuns:
% xdevf{j,k+1} is kth derivative of jth component).
%
% $Id: dde_sdmf_symdiff.m 170 2017-03-05 03:39:50Z jansieber $
%%
default={'maxorder',5,...
    'dev_append','_dev',...
    'deviation_name','h_devsmall',...
    'time_name','time',...
    'xxfunc_append','_func'};
options=dde_set_options(default,varargin,'pass_on');
n=size(xxs,1);
fs=fs(:);
taus=[0;taus(:)];
ntaus=length(taus);
%% define deviations
xnames=arrayfun(@(x)char(x),xxs(:,1),'uniformoutput',false);
xnamesrep=repmat(xnames,1,options.maxorder+1);
dt_inds=repmat((0:options.maxorder),[n,1]);
xdevnames=cellfun(@(x,y)[x,options.dev_append,'_',num2str(y)],...
    xnamesrep,num2cell(dt_inds),'uniformoutput',false);
xnamesrep=cellfun(@(x,y)[x,'_',num2str(y)],...
    xnamesrep,num2cell(dt_inds),'uniformoutput',false);
%% time argument t is argument of functions 
% xf and xdevf are functions of time, xdevf{j,k}(t) is k-1th derivative of jth
% component, xf{j}(t) is jth component
tname=options.time_name;
t=sym(tname);
xf=cellfun(@(x)symfun([x,'(',tname,')'],t),xnamesrep,'uniformoutput',false);
xdevf=cellfun(@(x)symfun([x,'(',tname,')'],t),xdevnames,'uniformoutput',false);
hname=options.deviation_name;
h=sym(hname);
if ~isempty(intersect(xnamesrep(:),xdevnames(:))) || ismember(hname,[xnamesrep(:);xdevnames(:)])
    error('symdiff:names',...
        'symdiff: name clash: choose different option for append or deviation_name\n%s\n%s\n%s\n',...
        evalc('celldisp(xnamesrep,''variables'')'),...
        evalc('celldisp(xdevnames,''deviations'')'),...
        evalc('celldisp(xnamesrep,''hname'')'));
end
%% substitute functions for arguments
% recursively insert values of delays into function xf{j}(t) and into
% right-hand side. Deviations are already included here.
ffunc=fs;
for k=1:ntaus
    for i=1:n
        for l=1:k-1
            taus(k)=subs(taus(k),xxs(i,l),...
                xf{i,1}(-taus(l))+h*xdevf{i,1}(-taus(l)));
        end
        ffunc=subs(ffunc,xxs(i,k),xf{i,1}(-taus(k))+h*xdevf{i,1}(-taus(k)));
    end
end
%% Differentiation wrt h
df=cell(1,options.maxorder);
for i=1:options.maxorder
    hd=repmat({h},1,i);
    df{i}=diff(ffunc,hd{:});
    df{i}=subs(df{i},h,0);
end
%% Set derivatives of base (equilibrium) to 0
% dfeq=df;
% for i=1:options.maxorder
%     for k=1:n
%         dfeq{i}=subs(dfeq{i},xf{k},xxs(k,1));
%     end
% end
%% replace kth derivative of xdevf{j,1} with xdevf{j,k+1}
% This is a hack, requiring a call to symengine, conversion of sym to
% character and manual interpretation (since expressions involving D(x)(t)
% are not substitutable by subs(...)
dfsubs=df;
for i=1:options.maxorder
    subseqs=cell(2,n);
    for k=1:n
        subseqs(:,k)={['D(',xdevnames{k,i},')=',xdevnames{k,i+1}];...
            ['D(',xnamesrep{k,i},')=',xnamesrep{k,i+1}]};
    end
    for l=1:options.maxorder
        dfsubs{l}=feval(symengine,'subs',char(dfsubs{l}),subseqs{:});
    end
end
%% Replace derivatives of x with zeros (in equilibrium)
% This is a hack as it substitutes the strings directly instead of the sym
% variables.
dfeq=dfsubs;
for l=1:options.maxorder
    for k=1:n
        dfeq{l}=subs(dfeq{l},xnamesrep{k,1},xnames{k});
        for i=2:options.maxorder+1
            dfeq{l}=subs(dfeq{l},xnamesrep{k,i},0);
        end
    end
end
end