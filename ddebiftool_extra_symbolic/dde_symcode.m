function str=dde_symcode(df,x,p,dx,dp,varargin)
%% convert symbolic expressions to matlab code using matlabFunction (uses temporary files!)
%
% output is string containing function body with name chooseable by
% optional input funcname.
%
% set optional input scalar to true to return every component as separate
% function (for sys_tau and its derivatives)
%
% $Id: dde_symcode.m 309 2018-10-28 19:02:42Z jansieber $
%%
%#ok<*AGROW>
default={'scalar',false,'funcname','sys','keeptemp',false,'outname','out',...
    'directional_derivative',true};
options=dde_set_options(default,varargin,'pass_on');
%% generate code with matlabFunction
% matlabFunction can only write to files, so generate temporary files where
% we save the function then read the file as string, concatenated 
nl=sprintf('\n');
str='';
if ~dde_isoctave()
    folder=tempname;
    gotfolder=mkdir(folder);
    if ~gotfolder
        error('dde_symcode:perm','dde_symcode: could not create temp folder %s',folder);
    end
end
%% break up output array into sequence of outputs 
% to enable scalar expansion after function call
vars={};
if ~options.directional_derivative
    vars{1}=[x(:).',p(:).'];
    vnames{1}=dde_names_from_sym(vars{1});
    for i=2:numel(df)
        dxd=horzcat(dx{1:i-1});
        dpd=horzcat(dp{1:i-1});
        vars{i}=[x(:).',p(:).',dxd(:).',dpd(:).'];
        vnames{i}=dde_names_from_sym(vars{i});
    end
else
    vars{1}=[x(:).',p(:).',dx{1}(:).',dp{1}(:).'];
    vnames{1}=dde_names_from_sym(vars{1});
end
nf=length(df{1});
if dde_isoctave()
    %% octave sym cannot handle setting output names
    options.out='out';
end
if ~options.scalar
    outnames=arrayfun(@(i)sprintf('%s_%d',options.outname,i),1:nf,'uniformoutput',false);
else
    outnames={options.outname};
    df=num2cell([df{:}]);
end
if dde_isoctave()
    outargs={'outputs',outnames};
else
    outargs={};
end
%%
if ~isempty(intersect([vnames{:}],outnames))
    error('dde_symcode:names',...
        'dde_symcode: name clash between output names and variables');
end
%% generate code in files tempname/funcname_ind_(order-1).m and re-read into string str
for i=1:numel(df)
    [ir,ic]=ind2sub(size(df),i);
    dfcell=num2cell(df{i});
    ind=1+(i-1)*(~options.directional_derivative);
    fname=sprintf('%s_%d_%d',options.funcname,ir,ic-1);
    if dde_isoctave()
        strnew=dde_octave_code(dfcell,fname,vars{ind});
    else
        % write code to temporary file
        filename=fullfile(folder,[fname,'.m']);
        w_orig=warning;
        warning('off','symbolic:generate:FunctionNotVerifiedToBeValid');
        matlabFunction(dfcell{:},'file',filename,'vars',vars{ind},outargs{:});
        warning(w_orig);
        % read code back into string str
        fid=fopen(filename,'r');
        strnew=fread(fid,inf);
        fclose(fid);
    end
    str=[str,nl,char(strnew(:)'),nl];
end
%% remove folder (unless optionally prevented)
if ~dde_isoctave && ~options.keeptemp
    rmdir(folder,'s')
end
end
