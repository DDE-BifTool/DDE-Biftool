%% Normal form calculations along Hopf bifurcation curves
% Using the normal form extension |nmfm| by Bram Wage.
%
% *Warning* Automatic bifurcation detection along branches is still in
% development. Its interface is likely to change, and detection is not yet
% reliable. With |br_bifdet|, normal form computations do not return error
% estimates. Thus, they should only be used with user-provided |sys_mfderi|
% for higher-order derivatives (in contrast to the default
% finite-difference approximation |mf_deriv|). See <demo1_simple.html> for
% an illustration, how one can perform a subset of the |br_bifdet|
% functionality with finite-difference approximations and error estimates.
%
% <html>
% $Id: demo1_normalforms.m 367 2019-07-14 21:56:46Z jansieber $
% </html>
%
% This demo requires to run <demo1_hopf.html> first.
%%
%#ok<*ASGLU,*NOPTS,*NASGU>
%
%% Load folder with normal form tools into matlab path and previous results
addpath('../../ddebiftool_extra_nmfm/',...
	'../../ddebiftool_utilities');
%% Higher-order derivatives provided by user
% Normal form calculations require higher-order derivatives. The user can
% provide them in the form of a function
%
%   function y=sys_mfderi(xx,par,v1,v2,...)
% 
% where |y=D^m_1f(xx,p)[v1,v2,...]|. The inputs |xx|, |v1|, |v2| are n x (ntau+1)
% arrays, |par| the system parameters. The output is the directional
% derivative of order |m|, where m is the the number of arguments |vk|
% (|k=1..m|).
%
% For this example, an explicit function for sys_mfderi is given, so we
% use it by adding the function to the function structure (not used below):
afuncs=set_funcs(funcs,'sys_mfderi',@sys_mfderi)
%% Automatic detection of codimension 2 bifurcations and criticality
% The |nmfm| extension provides |br_bifdet|, which determines Lyapunov
% cofficients, and codimension-two points along Hopf bifurcation curves.
%
% The output |branch_nmfm| has a |point| array of structures that contain
% additional fields |'nmfm'|, |'nvec'|, |'flag'|. The field |'flag'| may have the values
% 
% * |''| (empty): this point is not special,
% * |'hoho'|: this point is a double-Hopf point,
% * |'zeho'|: this point is a zero-Hopf point,
% * |'genh'|: this point is a generalized (degenerate/Bautin) Hopf point.
%
% The field |'nmfm'| is a structure. If the branch was a Hopf curve, it
% contains the field
%
% * |'L1'|: the first (3rd-order) Lyapunov coefficient at the Hopf point
% (|L1<0| implies that the Hopf bifurcation is supercritical, |L1>0| implies
% that the Hopf bifuration is subcritical).
% * for special points it contains other normal form coefficients,
% determining the sub-case of the bifurcation. See Kuznetsov (2004)
% 'Elements of Applied Bifurcation Theory') for notation and definitions.
%
[branch2_nmfm,br2_testfuncs,indbif_br2,types_br2]=LocateSpecialPoints(funcs,branch2);
%% Check Criticality
% The Lyapunov coefficients are always negative along the branch, implying
% that the bifurcation is supercritical. The are output as an array among
% the test functions for the generalized Hopf bifurcation.
L1_br2=br2_testfuncs.genh(1,:);
taus_br2=arrayfun(@(x)x.parameter(ind_taus),branch2_nmfm.point);
a21_br2=arrayfun(@(x)x.parameter(ind_a21),branch2_nmfm.point);
figure(17);
subplot(2,1,1);
plot(taus_br2,a21_br2,'.-');
grid on
xlabel('\tau_s');
ylabel('a_{21}');
title('(Repeated) two-parameter Hopf curve branch2');
subplot(2,1,2);
plot(taus_br2,L1_br2,'.-');
set(gca,'ylim',[-1,0]);
grid on
xlabel('\tau_s');
ylabel('L1');
title('First Lyapunov coefficient along branch2');
%% Special points
% The detection routine |br_bifdet| finds special points, listed by
% |br_getflags|. Each row index corresponds to a type of point. Indices can
% be converted to point types (strings) with |num2bif|.
%
% The call to |br_getflags| shows that |br_bifdet| has found two Hopf-Hopf
% points. Double-checking the data in the |'nmfm'| field shows that one of
% them is genuine, the other is spurious (a non-semisimple Hopf eigenvalue
% when a21=0). A typical feature of DDEs is that apparently
% higher-codimension degeneracies can occur due to the finite number of
% delays and equations.
%
% The points at which the Lyapunov coefficient L1 becomes singular
% (infinity) are special points: one of them is the non-semisimple Hopf
% eigenvalue point, the other is the Takens-Bogdanov point (where the Hopf
% frequency |omega| passes through 0). Again, this is not a generic
% Takens-Bogdanov point, but a Hopf-Pitchfork interaction (the system has
% reflection symmetry and the equilibrium is symmetric (0,0). Thus, the
% generic detection method is singular, reporting (correctly) that both
% normal form coefficients are 0.
special_br2_alltypes=br_getflags(branch2_nmfm)
knowntypes=num2bif();
biftype=num2bif(find(any(br_getflags(branch2_nmfm)>0,2)))
special_br2=special_br2_alltypes(bif2num(biftype),:)
for i=1:length(special_br2);
    fprintf('\nPoint %d normal form parameters:\n',special_br2(i));
    disp(branch2_nmfm.point(special_br2(i)).nmfm);
    fprintf('\nPoint %d Eigenvalues:\n',special_br2(i));
    disp(branch2_nmfm.point(special_br2(i)).stability.l0);
end
%% Perform bifurcation detection along second Hopf branch
% Along |branch3| we perform normal form computations and codimension-2
% detection, too.
[branch3_nmfm,br3_testfuncs,indbif_br3,types_br3]=LocateSpecialPoints(funcs,branch3);
%% Plot of Lyapunov coefficients for both Hopf branches
% Again, at |a21=0| we have a non-semisimple Hopf eigenvalue such that the
% Lyapunov coefficient tends to infinity. Otherwise, the Lyapunov
% coefficient is negative such that the Hopf bifurcation is supercritical.
L1_br3=br3_testfuncs.genh(1,:);
a21_br3=arrayfun(@(x)x.parameter(ind_a21),branch3_nmfm.point);
taus_br3=arrayfun(@(x)x.parameter(ind_taus),branch3_nmfm.point);
figure(18);
subplot(2,2,1);
plot(a21_br2,taus_br2,'.-',a21_br3,taus_br3,'.-');
a21lim=get(gca,'xlim');
tauslim=get(gca,'ylim');
grid on
xlabel('a_{21}');
ylabel('\tau_s');
title('(Repeated) two-parameter Hopf curves');
subplot(2,2,2);
plot(L1_br2,taus_br2,'.-');
set(gca,'ylim',tauslim,'xlim',[-0.6,0]);
grid on
ylabel('\tau_s');
xlabel('L1');
title('First Lyapunov coefficient along branch2');
subplot(2,2,3);
plot(a21_br3,L1_br3,'.-','color',[0,0.5,0]);
set(gca,'ylim',[-0.1,0],'xlim',a21lim);
grid on
xlabel('a_{21}');
ylabel('L1');
title('First Lyapunov coefficient along branch3');
%% Special points along branch3
% Only the Hopf-Hopf point already known from |branch2| and the spurious
% Hopf-Hopf point at |a21=0| are detected as codimension-2 points.
special_br3_alltypes=br_getflags(branch3_nmfm)
biftype=num2bif(find(any(br_getflags(branch3_nmfm)>0,2)))
special_br3=special_br3_alltypes(bif2num(biftype),:)
for i=1:length(special_br3)
    fprintf('\nPoint %d normal form parameters:\n',special_br3(i));
    disp(branch3_nmfm.point(special_br3(i)).nmfm);
    fprintf('\nPoint %d Eigenvalues:\n',special_br3(i));
    disp(branch3_nmfm.point(special_br3(i)).stability.l0);
end

%% Save and continue with periodic orbit continuation: <demo1_psol.html>
% See also <demo1_simple.html> for an illustration, how to use
% |ddebiftool_utilities| to perform a subset of the above computations
% without user-provided derivatives.
save('demo1_normalforms_results.mat');
