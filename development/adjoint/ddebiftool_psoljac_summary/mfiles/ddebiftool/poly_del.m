function dp=poly_del(m,c)
%% values of derivative of lagrange polynomials through 0:1/m:1 at c
% function dp=poly_del(m,c)
% INPUT:
%       m lagrange degree
%       c evaluation point(s) in [0,1]
% OUTPUT:
%       dp values of derivative of lagrange polynomials through 0:1/m:1 at c

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000
%
% $Id: poly_del.m 19 2014-04-11 14:15:36Z jan.sieber $
%
%% changed by js for vectorisation (exploited in psol_jac)
dp=zeros(m+1,length(c));
if m==1
    dp(1,:)=-1;
    dp(2,:)=1;
elseif m==2
    dp(1,:)=4*(c-0.75);
    dp(2,:)=-8*(c-0.5);
    dp(3,:)=4*(c-0.25);
elseif m==3
    dp(1,:)=-13.5*(c-0.4742165769367915).*(c-0.8591167563965418);
    dp(2,:)=40.5*(c-0.2615831876594899).*(c-0.8495279234516212);
    dp(3,:)=-40.5*(c-0.1504720765483788).*(c-0.7384168123405100);
    dp(4,:)=13.5*(c-0.1408832436034581).*(c-0.5257834230632086);
elseif m==4
    dp(1,:)=42.666666666666667*(c-0.3454915028125263).*(c-0.625).*(c-0.9045084971874737);
    dp(2,:)=-170.666666666666667*(c-0.16841363416054181147).*(c-0.61727336112119537025).*(c-0.90181300471826281827);
    dp(3,:)=256*(c-0.1047152924789526).*(c-0.5).*(c-0.8952847075210475);
    dp(4,:)=-170.666666666666667*(c-0.09818699528173718172).*(c-0.38272663887880462974).*(c-0.83158636583945818852);
    dp(5,:)=42.666666666666667*(c-0.09549150281252627).*(c-0.375).*(c-0.6545084971874737);
else
    t=0:1/m:1;
    for j=1:m+1
        for k=1:m+1
            if k~=j
                f=ones(1,length(c));
                for l=1:m+1
                    if l~=k && l~=j
                        f=f.*(c-t(l))/(t(j)-t(l));
                    end
                end
                dp(j,:)=dp(j,:)+f/(t(j)-t(k));
            end
        end
    end
end
% reshape to maintain compatibility with original
if length(c)==1
    dp=dp';
end
end
