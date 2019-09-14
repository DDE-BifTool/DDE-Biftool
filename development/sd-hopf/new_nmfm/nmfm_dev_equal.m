function iseq=nmfm_dev_equal(f1,f2)
%% Check if two history functions are (exactly) equal
%
% $Id$
%%
iseq=false;
[n1,nv1]=size(f1.v);
[n2,nv2]=size(f2.v);
if n2~=n1 || nv1~=nv2
    return
end
if any(f1.v(:)~=f2.v(:)) || any(f1.lambda~=f2.lambda) || any(f1.t~=f2.t)
    return
end
iseq=true;
end
