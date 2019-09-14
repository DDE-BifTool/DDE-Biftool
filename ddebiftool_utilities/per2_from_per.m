function [per2,suc]=per2_from_per(funcs,psolbranch,ind,varargin)
%% branch off at period doubling (branch psolbranch, point number ind)
%
% Wrapper for newer function DoublePsol
%
% $Id: per2_from_per.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
[per2,suc]=DoublePsol(funcs,psolbranch,ind,varargin);
end