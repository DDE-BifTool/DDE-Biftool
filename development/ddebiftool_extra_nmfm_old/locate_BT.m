%% correct point on fold or Hopf curve to Takens-Bogdanov point
%%
function [bt_point,success,btfuncs]=locate_BT(funcs,branch,ind)
%%
% [new_bifpoint,success]=locate_BT(funcs,branch,ind)
%
%% Input
%  funcs: right-hand side, delays, etc
%  branch: branch along which candidate for Takens-Bogdanov point was
%   detected
%  ind: index of branch.point to be used as starting point
%
%% Output
%  bt_point: point of kind 'BT' (with fields x, q0, q1)
%  success: flag indicating success of Newton iteration (non-zero=successful)
%
% $Id: locate_BT.m 115 2015-09-02 03:42:31Z jansieber $
%
%%
bifcand=p_tobt(funcs,branch.point(ind(1)));
%sys_deri=@(xx,par,nx,np,v)df_deriv(...
%    struct('sys_rhs',@(xx,par)sys_rhs_BT(xx,par,funcs)),xx,par,nx,np,v);
sys_deri=@(xx,par,nx,np,v)sys_deri_BT(xx,par,nx,np,v,funcs);
btfuncs=set_funcs(funcs,...
    'sys_rhs',@(xx,par)sys_rhs_BT(xx,par,funcs),...
    'sys_deri',sys_deri,...
    'sys_cond',@(p)sys_cond_BT(p,@extract_from_BT,funcs.sys_cond),'x_vectorized',false);
%btfuncs.sys_deri2=sys_deri2;
bifcand=struct('kind','stst','parameter',bifcand.parameter,...
    'x',[bifcand.x;bifcand.q0;bifcand.q1]);
cormethod=branch.method.point;
cormethod.extra_condition=1;
cormethod.minimal_accuracy=branch.method.bifurcation.secant_tolerance;
[corrected_point,success]=p_correc(btfuncs,bifcand,...
    branch.parameter.free,[],cormethod);
if ~success
    bt_point=[];
    return
end
corrected_point=nmfm_trim_point(corrected_point,branch.point(ind(1)));
bt_point=p_tobt(funcs,extract_from_BT(corrected_point,'solution'));
bt_point.q0=getfield(extract_from_BT(corrected_point,'q0'),'x');
bt_point.q1=getfield(extract_from_BT(corrected_point,'q1'),'x');
end
%% extract components of BT point
function result=extract_from_BT(btpoint,component)
dim=length(btpoint.x)/3;
type={'kind','solution','nullvector','q0','generalized nullvector','q1'};
switch component
    case type{1} % kind
        result='BT';
    case type{2} % solution
        result=btpoint;
        result.kind='stst';
        result.x=result.x(1:dim);
    case type(3:4) % eigenvector
        result=btpoint;
        result.kind='stst';
        result.x=result.x(dim+(1:dim));
    case type(5:6) % generalized eigenvector
        result=btpoint;
        result.kind='stst';
        result.x=result.x(2*dim+(1:dim));
    otherwise
        fprintf('known component types:\n');
        for k=1:length(type)
            fprintf('%s\n',type{k});
        end
        result=[];
end
end
