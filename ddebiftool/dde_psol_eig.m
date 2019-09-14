function stability=dde_psol_eig(funcs,pt,method,varargin)
%% compute stability information for stst 
% INPUT:
%   funcs problem functions
%	pt solution point
%	method method parameters 
% OUTPUT:
%	stability stability information
%
% experimental: method.sparse permits use of eigs
%
% $Id: dde_psol_eig.m 369 2019-08-27 00:07:02Z jansieber $
%
default={'check',true,'geteigenfuncs',false,'fill',[],'closest',[],...
    'eigmatrix','full','delay_accuracy',1e-8,'collocation_parameters',[],...
    'max_number_of_eigenvalues',20,'minimal_modulus',0};
[options,pass_on]=dde_set_options(default,[{'method',method},varargin],'pass_on','method');
pt=dde_psol_create(pt);
if pt.period<=0
    warning('p_stabil:period','p_stabil: period=%g,=0, no Floquet multipliers computed',pt.period);
    stability.mu=[];
    return
end
[J,res,tT,extmesh]=dde_coll_jac(funcs,pt,[],'c',options.collocation_parameters,...
    'Dtmat',funcs.lhs_matrix(size(pt.profile,1)),...
    'matrix',options.eigmatrix,...
    pass_on{:},'wrapJ',false,'output','profile'); %#ok<ASGLU>
dim=size(pt.profile,1);
em=repmat(extmesh,dim,1);
%% cut off negative delay if they are just numerical inaccuracy
d_ac=options.delay_accuracy;
if ~isempty(d_ac) && min(tT(:))>-d_ac && min(tT(:))<0 && max(em)>1
    [ir,ic,vJ]=find(J(:,em>=1));
    Jadd1=accumarray([ir(ic>dim),mod(ic(ic>dim)-1,dim)+1],vJ(ic>dim),[size(J,1),dim]);
    J=J(:,em<=1);
    J(:,em==1)=J(:,em==1)+Jadd1;
    em=em(em<=1);
    extmesh=extmesh(extmesh<=1);
    %tT=max(tT,0);
end
%% compute eigenvalues and eigenfunctions
[Margs,MPfun]=dde_psol_monodromy(J,options.closest,options.eigmatrix);
isef=options.geteigenfuncs;
[s,ef]=feval(['dde_psol_mult_app_',options.eigmatrix],Margs,...
    'geteigenfuncs',isef,'method',method,'closest',options.closest,...
    'max_number_of_eigenvalues',options.max_number_of_eigenvalues,pass_on{:});
assert(size(ef,2)==length(s)|~isef);
%% remove NaNs and Infs
isval=isfinite(s);
s=s(isval);
ef=ef(:,isval&isef);
%% sort, ensuring that positive imaginary parts come first
[~,im]=sort(imag(s),'descend');
s=s(im);
if isef
    ef=ef(:,im);
end
if isempty(options.closest)
    [~,is]=sort(abs(s),'descend');
else
    [~,is]=sort(abs(s-options.closest(1)),'ascend');
end
s=s(is);
if isef
    ef=ef(:,is);
end
%% remove unwanted eigenvalues   
sel= abs(s)>options.minimal_modulus & (1:length(s))'<=options.max_number_of_eigenvalues;
s=s(sel);
ef=ef(:,sel&isef);
neig=length(s);
err=NaN(1,neig);
if isef
    %% expand (if necesssary) eigenfunction to period interval + delay time
    ef=MPfun(ef);
    if options.check==1
        y=reshape(ef,dim,size(em,2),neig);
        exsel=extmesh<=extmesh(end)-1;
        exm1=extmesh(exsel);
        for i=1:neig
            y1=dde_coll_eva(y(:,:,i),extmesh,exm1+1,pt.degree);
            y0=y(:,exsel,i);
            err(i)=max(norm(J*ef(:,i),'inf'),norm(y0(:)*s(i)-y1(:),'inf'));
        end
    else
        err=NaN(1,neig);
    end
    %% convert eigenvectors to point structures
    ef=ef(em>=0&em<=1,:);
    efp=dde_point_from_x([ef;zeros(1,neig)+pt.period],pt,[]);
else
    efp=repmat(pt,1,0);
end
%% fill entries in s with 'fill'
if ~isempty(options.fill)
    s=cat(1,s,repmat(options.fill,options.max_number_of_eigenvalues-length(s),1));
end
%%
stability=struct('mu',s,'eigenfuncs',efp,'err',err);
end
% if ~funcs.tp_del
%     mu=mult_app(funcs,p.period,p.profile,p.mesh,p.degree,rho,max_number,col,p.parameter);
% else
%     mu=mult_app(funcs,p.period,p.profile,p.mesh,p.degree,rho,max_number,col,p.parameter,d_ac);
% end
% if ~isempty(mu)
%     [ydum,index_vector]=sort(abs(mu)); %#ok<ASGLU>
%     stability.mu=mu(index_vector(length(index_vector):-1:1));
% else
%     stability.mu=[];
% end
%end