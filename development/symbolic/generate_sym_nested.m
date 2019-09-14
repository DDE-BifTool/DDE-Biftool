%% Create right-hand side for nested_demo using the symbolic toolbox
% Demonstration for problems with state-dependent delays.
%
% <html>
% $Id$
% </html>
%
% The equation is 
% 
% $$x'(t)=-x(t+\theta_k)$$
%
% where $\theta_0(t)=0$ and
% 
% $$ \theta_{j+1}(t)=-p_1-x(t+\theta_j(t))-p_2 x(t+\theta_j(t))^2 $$
% 
% resulting in a nesting level of delays of order |k|. The parameter |p(1)|
% controls the delay at the Hopf bifurcation, |p(2)| controls the
% criticiality of the Hopf bifurcation (at least for |k=3| this is
% demonstrated in <nested_demo.html>).
%
%% Add paths
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');
%% The problem definition using symbolic variables
% The symbolic expression f is the right-hand side. The expression tau is
% the delay (ith component is delay number i).
ntau=3;  % nesting level of delays
nx=1;    % dimension of the problem
xx=sym('xx',[nx,ntau+1]);
par=sym('p',[1,2]);
f=-xx(1,ntau+1);
for i=1:ntau
    tau(i)=par(1)+xx(1,i)+par(2)*xx(1,i)^2;
end
%% Generation of right-hand side, delay function and their derivatives
% All functions are written to a single file with name set by optional
% argument |'filename'|. The optional argument |'sd_delay'| informs
% |dde_sym2funcs| that this is a state-dependent delay problem. For sd-DDEs
% not only the derivatives of f and tau are needed but also (for normal
% form computations the directional derivatives of the functions $f(x_t)$.
%
% The generated function |'sym_nested_tau_3'| (for |ntau=3|) can be passed
% on to |set_symfuncs|.
filename=sprintf('sym_nested_ntau_%d',ntau);
[funcstr,derivs]=dde_sym2funcs(f,xx,par,'sd_delay',tau,'filename',filename);
