function [points,ind]=dde_point_from_x(x,template_inp,free_par_ind)
%% insert variable values into point structure
%
% $Id: dde_point_from_x.m 308 2018-10-28 15:08:12Z jansieber $
%%
template=feval(['dde_',template_inp(1).kind,'_create'],template_inp(1));
template=dde_trim_point(template,template_inp(1));
ind=dde_ind_from_point(template,free_par_ind);
sx=size(x);
ptshape=[sx(2:end),1];
points=repmat(template,ptshape);
fnames=setdiff(fieldnames(ind),'parameter');
if isempty(x)
    return
end
npts=prod(ptshape);
x=reshape(x,size(x,1),npts);
for i=1:npts
    points(i).parameter(free_par_ind)=x(ind.parameter,i);
    for k=1:length(fnames)
        fd=ind.(fnames{k});
        if ~isstruct(fd)
            points(i).(fnames{k})=reshape(x(fd,i),size(template.(fnames{k})));
        else
            y=x(reshape(fd.re,[],1),i)+1i*x(reshape(fd.im,[],1),i);
            points(i).(fnames{k})=reshape(y,size(template.(fnames{k})));
        end
    end
end
end
