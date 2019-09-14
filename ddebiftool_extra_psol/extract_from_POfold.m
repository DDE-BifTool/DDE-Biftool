function result_array=extract_from_POfold(pfold_array,component,ip)
%% extract components from pfold solution branch or point array
%
% $Id: extract_from_POfold.m 369 2019-08-27 00:07:02Z jansieber $
%
%% check if input is branch rather than point array
if ~isfield(pfold_array,'kind') && isfield(pfold_array,'point')
    pfold_array=pfold_array.point;
end
%% extract named components
dim=ip.dim;
npar=ip.nuserpar;
type={'kind','solution','nullvector','solution_for_stability','xtau_ind'};
for i=1:length(pfold_array)
    pfold=pfold_array(i);
    switch component
        case type{1} %'kind'
            result_array='POfold';
            break
        case type{2} %'solution'
            result=pfold;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:npar);
        case type{3} %'nullvector'
            result=pfold;
            result.profile=result.profile(dim+1:end,:);
            result.parameter=result.parameter(ip.nullparind(:,2));
        case type{4} % including flag
            result=pfold;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:npar);
            result.flag='POfold';
        case type{5} % xtau_ind
            result_array=ip.ext_tau;
            break
        otherwise
            fprintf('known component types:\n');
            for k=1:length(type)
                fprintf('%s\n',type{k});
            end
            result_array=[];
            break;
    end
    result_array(i)=result; %#ok<AGROW>
end
end
