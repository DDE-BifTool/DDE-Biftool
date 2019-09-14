function [s,d,c]=p_dststb(point)
%% check stability and return dominant frequency and distance
% this function was apparently present in some original version of
% DDE-BifTool
%
% input: point (kind is stst or psol)
% output: 
% s stable (true) or not (false)
% c imaginary part of dominant eigenvalue
% d real part of dominant eigenvalue
%
% $Id: p_dststb.m 309 2018-10-28 19:02:42Z jansieber $
%

s=1;
c=0;
d=-Inf;
if ~isfield(point,'stability')
    error('p_dststb:stability', 'stability of point(s) not yet computed');
end
switch point.kind
    case 'stst'
        ev=point.stability.l0;
        if ~isempty(ev)
            s=sum(real(ev)>0)+1;
            [rdom,idom]=min(abs(real(ev))); %#ok<ASGLU>
            d=real(ev(idom));
            c=abs(imag(ev(idom)))>0;
        end
    case 'psol'
        ev=point.stability.mu;
        % remove trivial Floquet multiplier
        [edum,itriv]=min(abs(ev-1)); %#ok<ASGLU>
        ev(itriv)=[];
        if ~isempty(ev)
            s=sum(abs(ev)>1)+1;
            [rdom,idom]=min(abs(abs(ev)-1)); %#ok<ASGLU>
            d=abs(ev(idom))-1;
            if abs(imag(ev(idom)))>0
                c=1;
            elseif real(ev(idom))<0
                c=-1;
            else
                c=0;
            end
        end            
    otherwise
        error('p_dststb:kind', 'points of kind %s not supported',point.kind);
end
end
