function out=dde_postprocess(point,data)
%% postprocess a point after conversion from vector (if reqeuested by method)
% $Id: dde_postprocess.m 358 2019-07-04 22:46:43Z jansieber $
%%
if isfield(data.method,'postprocess') && ~isempty(data.method.postprocess)
    for i=numel(point):-1:1
        out(i)=feval(data.method.postprocess,point(i),data);
    end
    if numel(point)>0
        out=reshape(out,size(point));
    else
        out=point;
    end
else
    out=point;
end
end
