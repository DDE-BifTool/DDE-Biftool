function hcli=p_tohcli(point,varargin)
%% convert point to connecting orbit
% INPUT:
%     point a periodic solution near a homoclinic solution
%           alternatively an initial point in a hcli structure,
%           where a good starting guess for the profile and steady
%           states are available
%     named (but mandatory) 'funcs': problem functions
% OUTPUT:
%     hcli a starting value to compute the exact homoclinic or
%     heteroclinic solution  
%
% (c) DDE-BIFTOOL v. 2.02, 16/6/2002
%
% $Id: p_tohcli.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
%% for backward compatibility re-organize arguments
args=varargin;
if ~isfield(point,'kind')
    funcs=point;
    point=args{1};
    args=[args(2:end),{'funcs',funcs}];
end
default={'funcs',[]};
options=dde_set_options(default,args,'pass_on'); 
sys_tau=options.funcs.sys_tau;
sys_deri=options.funcs.sys_deri;

if mod(length(point.mesh),point.degree)~=1
    err=[length(point.mesh) point.degree];
    error('P_TOHCLI: psol does not contain L=%d intervals of m=% points!',...
        err(1),err(2));
end

hcli.kind='hcli';
hcli.parameter=point.parameter;
hcli.mesh=point.mesh;
hcli.degree=point.degree;

switch point.kind
    case 'psol'
        ntst=size(point.profile,2);
        test=NaN(1,ntst-1);
        for i=1:ntst-1
            test(i)=norm(point.profile(:,i)-point.profile(:,i+1));
        end
        [minval, pos]=min(abs(test)); %#ok<ASGLU>
        stst.kind='stst';
        stst.parameter=hcli.parameter;
        stst.x=point.profile(:,pos);
        x_profile=NaN(1,ntst);
        for i=1:size(point.profile,2)
            x_profile(1,i)=norm(point.profile(:,i)-stst.x);
        end
        
        [peak, peak_pos]=max(x_profile); %#ok<ASGLU>
        [hole, hole_pos]=min(x_profile); %#ok<ASGLU>
        left_part=point.profile(:,1:peak_pos);
        right_part=point.profile(:,peak_pos+1:end);
        hole_begin=hole_pos-mod(hole_pos,point.degree)+1;
        hole_end=hole_begin+point.degree;
        
        if hole_pos<peak_pos
            right_part=[right_part left_part(:,2:hole_begin)];
            left_part=left_part(:,hole_end:end);
            hcli.mesh=[hcli.mesh(hole_end:end) ...
                (hcli.mesh(2:hole_begin)+1)];
        else
            left_part=[right_part(:,hole_end-peak_pos:end-1) left_part];
            right_part=right_part(:,1:hole_begin-peak_pos);
            hcli.mesh= ...
                [(hcli.mesh(hole_end:end-1)-1)...
                hcli.mesh(1:hole_begin)];
        end
        
        nb_of_points=length(hcli.mesh);
        rest=mod(nb_of_points,point.degree);
        hcli.profile=[left_part right_part];
        
        if rest>1
            hcli.profile=point.profile(:,1+floor((rest-1)/2):end-ceil((rest-1)/2));
            hcli.mesh=hcli.mesh(1+floor((rest-1)/2):end-ceil((rest-1)/2));
        end
        if rest==0
            rest=point.degree;
            hcli.profile=point.profile(:,1+floor((rest-1)/2):end-ceil((rest-1)/2));
            hcli.mesh=hcli.mesh(1+floor((rest-1)/2):end-ceil((rest-1)/2));
        end
        hcli.mesh=hcli.mesh-hcli.mesh(1);
        hcli.period=point.period*hcli.mesh(end);
        hcli.mesh=hcli.mesh/hcli.mesh(end);
        hcli.x1=stst.x;
        hcli.x2=stst.x;
        stst1=stst;
        stst2=stst;
    case 'hcli'
        hcli=point;
        stst1.kind='stst';
        stst1.parameter=point.parameter;
        stst1.x=point.x1;
        stst2=stst1;
        stst2.x=point.x2;
    otherwise
        error(['P_TOHCLI: not a valid conversion for other than psol' ...
            ' or hcli type points']);
end

m=df_mthod('stst','cheb');
stst1.stability=p_stabil(options.funcs,stst1,m.stability);

if isempty(stst1.stability.l1) || max(real(stst1.stability.l1))<0
    error('P_TOHCLI: no unstable eigenmodes found');
end
lambda_v=stst1.stability.l1(:);
v_sel=real(lambda_v)>0;
lambda_v=lambda_v(v_sel);

hcli.lambda_v=lambda_v;
hcli.v=stst1.stability.v(:,v_sel);

stst2.stability=p_stabil(funcs,stst2,m.stability);

lambda_w=stst2.stability.l1(:);
w_sel=real(lambda_w)>0;
lambda_w=lambda_w(w_sel);
hcli.lambda_w=lambda_w;
hcli.w=stst2.stability.w(:,w_sel);

hcli.alpha=hcli.v\(hcli.profile(:,1)-hcli.x1);
hcli.epsilon=norm(hcli.alpha);
hcli.alpha=hcli.alpha/hcli.epsilon;
end
