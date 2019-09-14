function [pout,conversion]=dde_coll_convert(point,conversion)
%% combine coll time profiles and parameters for extended system to single profiles 
% (used for psol) and its extensions
%
% $Id: dde_coll_convert.m 339 2019-05-09 19:47:01Z jansieber $
%%
if isempty(point)
    pout=point;
    return
end
if nargin>1 && ~isempty(conversion)
    pout=feval(['dde_',conversion.kind,'_create'],'point',point);
    fnames=setdiff(fieldnames(pout),fieldnames(point));
    for i=1:length(fnames)
        pout.(fnames{i})=conversion.copy.(fnames{i});
    end
    pout=add_fields(pout,conversion,point,'profile',1);
    pout=add_fields(pout,conversion,point,'parameter',2);
else
    conversion.kind=point.kind;
    [dum,dum,convtarget]=feval(['dde_',conversion.kind,'_create'],'point',point); %#ok<ASGLU>
    pout=dde_coll_create('point',point);
    pout.kind=convtarget;
    [pout,conversion]=reduce_fields(pout,conversion,point,'profile',1);
    [pout,conversion]=reduce_fields(pout,conversion,point,'parameter',2);
    conversion.copy=point;
end
end
%%
function [pout,conversion]=reduce_fields(pout,conversion,point,name,dim)
[~,varout]=feval(['dde_',point.kind,'_create']);
if iscell(varout)
    varfields=varout{1};
else
    varfields=varout;
end
fn=fieldnames(varfields);
fn=fn(strncmp(fn,[name,'_'],length(name)+1));
ind0=size(point.(name),dim);
ind=ind0;
for i=1:length(fn)
    len=size(point.(fn{i}),dim);
    conversion.(fn{i})=ind+(1:len);
    ind=ind+len;
end
conversion.extra_ind.(name)=ind0+1:ind;
for i=length(fn):-1:1
    ar=pout.(name);
    if dim==1
        ar(conversion.(fn{i}),:)=point.(fn{i});
    elseif dim==2
        ar(:,conversion.(fn{i}))=point.(fn{i});
    end
    pout.(name)=ar;
end
end
%%
function pout=add_fields(pout,conversion,point,name,dim)
fields=fieldnames(conversion);
fn=fields(strncmp(fields,[name,'_'],length(name)+1));
for i=1:length(fn)
    ar=point.(name);
    if dim==1
        pout.(fn{i})=ar(conversion.(fn{i}),:);
    elseif dim==2
        pout.(fn{i})=ar(:,conversion.(fn{i}));
    end
end
ar=point.(name);
if dim==1
    ar(conversion.extra_ind.(name),:)=[];
elseif dim==2
    ar(conversion.extra_ind.(name))=[];
end
pout.(name)=ar;
end
