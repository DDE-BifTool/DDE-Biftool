function method=dde_psolfold_method(varargin)
%% define method parameters for psolfold
% $Id: dde_psolfold_method.m 362 2019-07-14 15:49:40Z jansieber $
%% 
if isstruct(varargin{1})
    method=varargin{1};
    args=varargin(2:end);
else
    method=df_mthod('psol');
    args=varargin;
end
method.point.hdev=0;
method=replace_method_pars(method,args{:});
method.point.preprocess='dde_psolfold_preprocess';
method.point.postprocess='dde_psolfold_postprocess';
end
