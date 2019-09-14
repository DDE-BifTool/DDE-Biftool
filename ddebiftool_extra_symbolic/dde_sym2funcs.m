function [funcstr,derivatives]=dde_sym2funcs(fs,xxs,ps,varargin)
%% create matlab code for r.h.s (and delays if needed) from symbolic expressions
%
%   function [funcstr,derivatives]=dde_sym2funcs(fs,xxs,ps,...) 
%
% It converts the input fs (an array of symbolc expressions) into a
% right-hand side function and its derivatives that can be used for
% DDE-Biftool. The wrapper |set_symfuncs| around |set_funcs| can be used to
% create a structure containing right-hand sides and derivatives.
%
%% Inputs:
%
% * |fs|:  n x 1 array of symbolic expressions, defining the right-hand side
% * |xxs|: n x (ntau+1) array of symbols for states, |xxs(j,k)| is the
% $x_j(t-\tau_{k-1})$, and |xxs(j,1)|=$x_j(t)$.
% * |par|: 1 x np (or np x 1) array of symbols for parameters, |par(j)| is
% the |j|the parameter.
%
%% Common optional name-value input pairs
%
% * |'sd_delay'| (default |sym([]))|): ntau x 1 (or 1 x ntau) array tau of
% symbolic expressions depending on symbols in |xxs| and |par|. |tau(k)| is
% delay number |k| ($\tau_{k}(x,p)$) an dmay depend on |xxs(j,l)| for
% |l=1..k|.
% * |'write'| (default |true|): write output to file?
% * |'filename'| (default |'sys'|): results are written to function file
% with this name.
% * |'maxorder'| (default 5): maximal order of derivatives to be computed.
%
%% Outputs (often not needed)
%
% * |fstr|: strong containing the code for writing the function.
% * |derivatives|: symbolic expressions for all dervatives (a structure,
% containing the fields |df|, |xx|, |parameter|, |dx| and |dp|. |df|
% contains the expressions, |xx|, |parameter|, |dx| and |dp| the symbols
% used in these expressions. 
%
% The routine will create a vectorized matlab
% function |y=f(xx,p)| from the expressions in |fs|, and its directional
% derivatives up to a maximum order (default 5). For state-dependent delays
% (optional argument |sd_delay| non-empty), it also computes the values and
% directional derivatives of all delays.
%% Warning
% The file will write to temporary files, since matlabFunction can only
% write to files.
%
% $Id: dde_sym2funcs.m 320 2019-02-01 00:38:43Z jansieber $
%%
default={'svn_id','','write',true,'filename','sys','sd_delay',sym([]),'directional_derivative',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% use directional derivative if delay is state dependent
% because multilinear derivatives are not yet implemented
if ~isempty(options.sd_delay)
    options.directional_derivative=true;
end
pass_on=[{'directional_derivative',options.directional_derivative},pass_on];
%% convert possibly cell arrays safely into sym arrays
fs=dde_sym_from_cell(fs);
xxs=dde_sym_from_cell(xxs);
ps=dde_sym_from_cell(ps);
options.sd_delay=dde_sym_from_cell(options.sd_delay);
fs=fs(:);
%% extract function name from (optional) filename
[pth,funcname,ext]=fileparts(options.filename);
[fstr{1},df{1},v{1},q{1}]=dde_symdericode(fs,xxs,ps,[funcname,'_rhs'],pass_on{:});
maxorder=length(df{1})-1;
if ~isempty(options.sd_delay)
    tp_del=true;
    tau=options.sd_delay;
    ntau=length(tau);
    [fstr{2},df{2},v{2},q{2}]=dde_symdericode(tau,xxs,ps,[funcname,'_tau'],'scalar',true,pass_on{:});
    %[fstr{3},df{3},v{3}]=dde_sdmf_symdericode(fs,tau,xxs,ps,[funcname,'_combined'],pass_on{:});
    %mf_dxlength=numel(v{3});
else
    ntau=size(xxs,2)-1;
    tp_del=false;
    %mf_dxlength=0;
end
str=[fstr{:}];
derivatives=struct('df',df,'xx',xxs,'parameter',ps,'dx',v,'dp',q);
%% octave's output ends with 'end', matlab's does not
if dde_isoctave()
    function_end='end';
else
    function_end='';
end
%% create full function (still string)
nl=sprintf('\n');
header=sprintf('function varargout=%s(action,varargin)',funcname);
comment=[...
    '%% Automatically generated with matlabFunction',nl,...
    '% ',options.svn_id,nl,...
    '%#ok<*DEFNU,*INUSD,*INUSL>',nl];
body=[...
    'switch action',nl,...
    '  case ''ntau''',nl,...
    '   varargout{1}=',num2str(ntau),';',nl,...
    '   return',nl,...
    '  case ''tp_del''',nl,...
    '   varargout{1}=',num2str(tp_del),';',nl,...
    '   return',nl,...    
    '  case ''maxorder''',nl,...
    '   varargout{1}=',num2str(maxorder),';',nl,...
    '   return',nl,... 
    '  case ''directional_derivative''',nl,...
    '   varargout{1}=',num2str(options.directional_derivative),';',nl,...
    '   return',nl,... 
    'end',nl,...
    'ind=varargin{1};',nl,...
    'order=varargin{2};',nl,...
    'nout=varargin{3};',nl,...
    'f=str2func(sprintf(''',funcname,'_%s_%d_%d'',action,ind,order));',nl,...
    'varargout=cell(nout,1);',nl,...
    '[varargout{:}]=f(varargin{4:end});',nl,...
    function_end,nl,...
    nl];
funcstr=[header,nl,comment,nl,body,nl,str];
%% write to file if requested
if options.write
    if isempty(ext)
        ext='.m';
    end
    filename=fullfile(pth,[funcname,ext]);
    fid=fopen(filename,'w');
    if fid<=0
        error('dde_sym2funcs:perm','dde_sym2funcs: could not create function file %s',filename);
    end
    fwrite(fid,funcstr);
    fclose(fid);
end
