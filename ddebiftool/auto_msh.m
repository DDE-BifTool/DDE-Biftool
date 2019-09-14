function tmnew=auto_msh(ups,tmold,nnew)
%% create new mesh on [0,1] equidistributing the error
% function tmnew=auto_msh(ups,tmold,dtmold,nnew)
% INPUT:
%	ups solution profile
%	tmold old mesh
%	nnew new number of mesh points
% OUTPUT:
%	tmnew new mesh
% COMMENT: 
%	this function is a matlab translation of the AUTO
%	fortran function NEWMSH restricted to the case of periodic 
%	boundary conditions

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000
%
% $Id: auto_msh.m 296 2018-09-24 21:36:56Z jansieber $
%
%% put the values of the monotonely increasing error in eqf:
dtmold=tmold(2:end)-tmold(1:end-1);
eqf=auto_eqd(dtmold,ups);

%% uniformly divide the range of eqf:
nnewp1=nnew+1;
uneq=linspace(0,eqf(end),nnew+1);
%% find how indices have to be shifted
ial=order_indices(eqf,uneq);

%% generate the new mesh in tmnew
x=(uneq-eqf(ial))./(eqf(ial+1)-eqf(ial));
tmnew=(1-x).*tmold(ial)+x.*tmold(ial+1);
tmnew(nnewp1)=1;
end
%% find interval indices in mesh tm for points in tm1
function itm1=order_indices(tm,tm1)
% function itm1=order_indices(tm,tm1)
% INPUT:
%	tm ascending array 
%	tm1 ascending array
% OUTPUT:
%	itm1 itm1(k) gives index of tm1(k) in tm-interval
%
%%
n=length(tm);
n1=length(tm1);
itm1=NaN(1,n1);
[tdum,itx]=sort([tm,tm1]); %#ok<ASGLU>
jt=0;
jx=1;
for i=1:length(itx)
    if itx(i)<=n
        jt=jt+1;
    else
        itm1(jx)=jt;
        jx=jx+1;
    end
end
if tm(end)==tm1(end)
    itm1(end)=itm1(end)-1;
end
end
