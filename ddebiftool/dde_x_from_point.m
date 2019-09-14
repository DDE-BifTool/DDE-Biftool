function [x,ind]=dde_x_from_point(points,free_par_ind)
%% extract variable values from point structure
%
% $Id: dde_x_from_point.m 374 2019-09-14 14:02:58Z jansieber $
%%
x=[];
ind=repmat(struct(),1,0);
if isempty(points)
    return
end
[ind,len]=dde_ind_from_point(points(1),free_par_ind);
x=NaN(len,numel(points));
fnames=fieldnames(ind);
for i=1:length(points)
    points(i).parameter=points(i).parameter(free_par_ind);
    for k=1:length(fnames)
        fd=ind.(fnames{k});
        if ~isstruct(fd)
            x(reshape(fd,[],1),i)=reshape(points(i).(fnames{k}),[],1);
        else
            x(reshape(fd.re,[],1),i)=reshape(real(points(i).(fnames{k})),[],1);
            x(reshape(fd.im,[],1),i)=reshape(imag(points(i).(fnames{k})),[],1);
        end
    end 
end
x=reshape(x,[len,size(points)]);
end
