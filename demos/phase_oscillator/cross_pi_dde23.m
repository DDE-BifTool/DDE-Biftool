%% set event for dde23
% generates event whenever y crosses |pi|
function [value,isterm,dir]=cross_pi_dde23(y)
value=sin(y);
isterm=false;
dir=1;
end
