function funcstr=dde_sdmf_symcode(df,x,p,dx,varargin)
%% convert symbolic expressions to matlab code using matlabFunction (uses temporary files!)
%
% output is string containing function body with name chooseable by
% optional input funcname.
%
%% Inputs
%
% * df: 1 x maxorder cell array of n x 1 sym arrays for the right-hand side
% df{j} is the jth order functional derivative of f(x+h*dx) wrt to h.
% x is an equilibrium, so a constant, but dx is a cell array of symfuns:
% dx{j,k+1} is kth derivative of jth component.
% * x: n x nc array of syms (nc should be 1 as this works only for
% equilibria)
% * p: 1 x npar vector of sym parameters (passed on as arguments to
% generated code)
% * dx  cell array of symfuns: dx{j,k+1} is kth derivative of jth
% component.
%
% $Id: dde_sdmf_symcode.m 168 2017-03-03 22:04:32Z jansieber $
%%
%#ok<*AGROW>
default={'funcname','sys_dirderi','keeptemp',false,'svn_id','','outname','out'};
options=dde_set_options(default,varargin,'pass_on');
%% generate code with matlabFunction
% matlabFunction can only write to files, so generate temporary files where
% we save the function then read the file as string, concatenated 
nl=sprintf('\n');
str='';
folder=tempname;
gotfolder=mkdir(folder);
if ~gotfolder
    error('dde_sdmf_symcode:perm','dde_sdmf_symcode: could not create temp folder %s',folder);
end
%% break up output array into sequence of outputs 
% to enable scalar expansion after function call
vars=[x(:).',p(:).',dx(:).'];
vnames=arrayfun(@(x)char(x),vars,'uniformoutput',false);
nf=length(df{1});
outnames=arrayfun(@(i)sprintf('%s_%d',options.outname,i),1:nf,'uniformoutput',false);
if ~isempty(intersect(vnames,outnames))
    error('dde_sdmf_symcode:names',...
        'dde_sdmf_symcode: name clash between output names and variables');
end
%% generate code in files tempname/funcname_ind_(order-1).m and re-read into string str
for i=1:numel(df)
    [ir,ic]=ind2sub(size(df),i);
    % write code to temporary file
    fname=fullfile(folder,sprintf('%s_%d_%d.m',options.funcname,ir,ic-1));
    dfcell=num2cell(df{i});
    matlabFunction(dfcell{:},'File',fname,'Vars',vars,'Outputs',outnames);
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
%% create full function (still string)
header=sprintf('function varargout=%s(ind,order,nout,varargin)\n',options.funcname);
comment=[...
    '%% Automatically generated with matlabFunction',nl,...
    '% ',options.svn_id,nl,...
    '%#ok<*DEFNU,*INUSD,*INUSL>',nl];
body=[sprintf(...
    'f=str2func(sprintf(''%s_',options.funcname),'%d_%d'',ind,order));',nl,...
    'varargout=cell(nout,1);',nl,...
    '[varargout{:}]=f(varargin{:});',nl,...
    nl];
funcstr=[header,comment,body,str];
end
