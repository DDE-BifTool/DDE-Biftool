function str=dde_symcode(df,x,p,dx,dp,varargin)
%% convert symbolic expressions to matlab code using matlabFunction (uses temporary files!)
%
% output is string containing function body with name chooseable by
% optional input funcname.
%
% set optional input scalar to true to return every component as separate
% function (for sys_tau and its derivatives)
%
% $Id: dde_symcode.m 170 2017-03-05 03:39:50Z jansieber $
%%
%#ok<*AGROW>
default={'scalar',false,'funcname','sys','keeptemp',false,'outname','out'};
options=dde_set_options(default,varargin,'pass_on');
%% generate code with matlabFunction
% matlabFunction can only write to files, so generate temporary files where
% we save the function then read the file as string, concatenated 
nl=sprintf('\n');
str='';
folder=tempname;
gotfolder=mkdir(folder);
if ~gotfolder
    error('dde_symcode:perm','dde_symcode: could not create temp folder %s',folder);
end
%% break up output array into sequence of outputs 
% to enable scalar expansion after function call
vars=[x(:).',p(:).',dx(:).',dp(:).'];
vnames=arrayfun(@(x)char(x),vars,'uniformoutput',false);
nf=length(df{1});
if ~options.scalar
    outnames=arrayfun(@(i)sprintf('%s_%d',options.outname,i),1:nf,'uniformoutput',false);
else
    outnames={options.outname};
    df=num2cell([df{:}]);
end
if ~isempty(intersect(vnames,outnames))
    error('dde_symcode:names',...
        'dde_symcode: name clash between output names and variables');
end
%% generate code in files tempname/funcname_ind_(order-1).m and re-read into string str
for i=1:numel(df)
    [ir,ic]=ind2sub(size(df),i);
    % write code to temporary file
    fname=fullfile(folder,sprintf('%s_%d_%d.m',options.funcname,ir,ic-1));
    dfcell=num2cell(df{i});
    w_orig=warning;
    warning('off','symbolic:generate:FunctionNotVerifiedToBeValid');
    matlabFunction(dfcell{:},'File',fname,'Vars',vars,'Outputs',outnames);
    warning(w_orig);
    % read code back into string str
    fid=fopen(fname,'r');
    strnew=fread(fid,inf);
    fclose(fid);
    str=[str,nl,char(strnew(:)'),nl];
end
%% remove folder (unless optionally prevented)
if ~options.keeptemp
    rmdir(folder,'s')
end
end
