function [C1funcs,C1branchout,suc] = C1branch_from_C2point(funcs,C2point,freepars,varargin)
%% Initialize branch for continuing the Limit point of cycles curve
% emanating from the the generalized-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% $Id: C1branch_from_C2point.m 314 2019-01-24 14:28:23Z mmbosschaert $
%%
default={'step',1e-2,'codim2','','codim1','','correc',true,'print',1,...
    'c1args',{},'predictor',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isempty(options.codim1)
    error('C1branch_from_C2point:args',...
        'C1branch_from_C2point: provide argument codim1');
end
if isempty(options.codim2)
    error('C1branch_from_C2point:args',...
        'C1branch_from_C2point: provide argument codim2');
end
if ~strcmp(C2point.kind,options.codim2)
    p_toc2=str2func(['p_to',options.codim2]);
    nmfm_c2=str2func(['nmfm_',options.codim2]);
    C2point=p_toc2(C2point);
    C2point=nmfm_c2(funcs,C2point,'free_pars',freepars);
end
%% approximate first (at least) two points
if length(options.step)==1
    options.step(2)=options.step(1)*2;
end
c1approx=str2func(['nmfm_',options.codim1,'_from_',options.codim2,'_init']);
[br,augmented]=c1approx(funcs,C2point,options.step,freepars,pass_on{:});
C1funcs=funcs;
C1branch=br;
if (islogical(augmented) && ~augmented) || iscell(augmented)
    SetupC1funcs=str2func(['Setup',options.codim1]);
    for k=length(br):-1:1
        if iscell(augmented) && ~isempty(augmented)
            optargs=augmented{k};
        else
            optargs={};
        end
        for i=1:length(options.step)
            [C1funcs,tmpbranch]=SetupC1funcs(funcs,br(k),i,pass_on{:},...
                'contpar',freepars,'correc',false,optargs{:});
            if i==1
                C1branch(k)=tmpbranch;
            else
                C1branch(k).point(i)=tmpbranch.point(1);
            end
        end
    end
else
    C1branch=br;
    for k=length(br):-1:1
        C1branch(k)=replace_branch_pars(br(k),freepars,pass_on);
    end
end
suc=false(length(options.step),length(C1branch));
for k=length(C1branch):-1:1
    if isfield(br(k),'tangent')
        stepcond=br(k).tangent(C1branch(k).point);
    else
        stepcond=p_axpy(-1,C1branch(k).point(1),C1branch(k).point(2));
        stepdist=p_norm(stepcond);
        stepcond=p_axpy(1/stepdist,stepcond,[]);
    end
    if ~options.predictor
        for i=1:length(options.step)
            refpoint=C1branch(k).point(max(i-1,1));
            [C1branch(k).point(i),suc(i,k)]=p_correc(C1funcs,C1branch(k).point(i),...
                C1branch(k).parameter.free,stepcond,C1branch(k).method.point,0,...
                refpoint);
            if options.print>0
                fprintf([options.codim1,' from ',options.codim2,...
                    ': branch %d of % d correction of point %d, success=%d\n'],...
                    k,length(C1branch),i,suc(i,k));
            end
        end
    end
    if isfield(C1branch(k),'tangent')
        C1branchout(k)=rmfield(C1branch(k),'tangent');
    else
        C1branchout(k)=C1branch(k);
    end
end
end
