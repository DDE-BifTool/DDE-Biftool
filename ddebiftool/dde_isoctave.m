function isoct=dde_isoctave()
%% check if octave is running instead of Matlab
%
% $Id: dde_isoctave.m 296 2018-09-24 21:36:56Z jansieber $
%
if exist('OCTAVE_VERSION','builtin')
    isoct=true;
else
    isoct=false;
end
end
