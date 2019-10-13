%% DDE-Biftool - Demonstration of linear stability analysis
%
% <html>
% $Id: linear_stablity_demo.m 368 2019-07-15 07:53:59Z jansieber $
% </html>
%
% This demo illustrates common ways to perform linear stability analysis of
% equilibria and periodic orbits and potential problems that are typical
% for DDEs. Numerical problems arise whenever the delay is "large" in some
% sense.
%% Add paths and clear variables
clear
addpath([pwd(),'/../../ddebiftool'],...
    [pwd(),'/../../ddebiftool_extra_psol'],...
    [pwd(),'/../../ddebiftool_extra_nmfm/'],...
    [pwd(),'/../../ddebiftool_utilities']);
format compact
format short g
%% Standard usage after computation of branches
% Consider the scalar DDE
% 
% $$x'(t)=-x(t-tau)+x(t)^3-x(t)^5$$
%
% This DDE has a unique equilibrium at $x=0$, which loses its stabilit in a
% subcritical Hopf bifurcation at $\tau=\pi/2$. The family of periodic
% orbits branching off is initially unstable but then becomes stable for
% larger amplitude.
%% Definition of right-hand side and computation of trivial equilibrium branch
funcs=set_funcs('sys_rhs',@(x,p)-x(1,2,:)+x(1,1,:).^3-x(1,1,:).^5,...
    'sys_tau',@()1,'x_vectorized',true);
[eqbr,suc]=SetupStst(funcs,'x',0,'parameter',0,'contpar',1,...
    'step',0.2,'max_bound',[1,pi])
figure(1);clf;
eqbr=br_contn(funcs,eqbr,20);
%% Stability of branches of equilibria with GetStability
% The utility function |GetStability| performs linear stability analysis on
% every point on the branch. It returns as its first output the number of
% unstable eigenvalues. This number takes into account the type of point
% (for example, it will discard the pair of purely imaginary eigenvalues
% for points on branches of Hopf bifurcations). The fourth output is the
% array of points, with stability information attached (second and third
% output can be ignored for equilibria). An alternative to |GetStability|
% is |br_stabl|, which only computes the stability information without
% counting the number of unstable eigenvalues. Both functions do not
% recompute stability information if it is already present, unless
% requested.
[nunst_eq,~,~,eqbr.point]=GetStability(eqbr,'funcs',funcs)
%% Point stability
% Inside |GetStability| or |br_stabl| the function |p_stabil| is called to
% populate the field |stability| of the |point| structure. The |stability|
% field's entries depend on the method used for stability computation.
% Parameters are listed in |method.stability|, which is passed on as a
% third argument of p_stabil. Parameters can be varied as optional
% arguments. Of particular interest are |minimal_real_part|, |maxsize| (to
% increase the maximal size of the discretized matrix),
% |max_number_of_eigenvalues|.
pt=eqbr.point(end);
stability=p_stabil(funcs,pt,eqbr.method.stability,...
    'max_number_of_eigenvalues',Inf,'minimal_real_part',-2,'maxsize',1000);
fprintf('Stability field\n');
disp(stability)
figure(1);clf;
p_splot(stability);
%% Lower-level functions
% Depending on the parameter |discretization| in |method.stability| a
% lower-level linear stability routine |dde_stst_eig_yyy| is called, where
% |yyy| can be either 
%
% * |cheb| (the default), based on Breda et al's traceDDE discretization
% using Chebyshev polynomials and the infinitesimal generator,
% * |bdf|, the original DDE-Biftool stability method by Engelborghs et
% al, based on the linear multistep (LMS) BDF method for a single time step of
% the linear semigroup,
% * |mxo|, the modification provided with a different heuristics by
% Verheyden et al (also a LMS method).
%
% These functions have a |n x n x (d+1)| array |A| of coefficient matrices
% as their first input and a |1xd| array |tau| as their second input. All
% other inputs are optional, one optional argument can be |'method'| with a
% method stability structure as created by |df_mthod| (and present by
% default in |eqbr.method|.
%
% The outputs are always structures containing the fields |l0| and |l1|,
% containing the eigenvalues. These fields may have slightly different
% entries for |dde_stst_eig_bdf| and |dde_stst_eig_mxo|, but are identical
% except for transposition for |dde_stst_eig_cheb|. All other fields depend
% on the method. The method |dde_stst_eig_cheb| also provides left and
% right eigenvectors in the fields |w| and |v|, and an error estimate for
% the eigenvalue |err|, which equals the norm of |Delta(l0(i))*v(:,i)|. The
% field |discarded|, if not empty, contains eigenvalues of the discretized
% problem for which this residual was larger than
% |root_accuracy*discard_accuracy_factor| (|1e-6*1e5=0.1| is the default).
% A non-empty field |discarded| points to numerical problems that require
% increasing the permitted |maxsize| or increasing |minimal_real_part|
% (pushing it closer to |0|).
%
% The above result can be obtained with continuation and creation of
% right-hand side or point stuctures for |x'=-x(t-tau)| and |tau=pi| using
stability=dde_stst_eig_cheb(cat(3,0,-1),pi,...
    'max_number_of_eigenvalues',Inf,'minimal_real_part',-2,'maxsize',1000);
figure(1);clf;
p_splot(stability);
hold on
pli=plot(stability.discarded.l0,'bx','DisplayName','inaccurate');
legend(pli)
%% Problems with accuracy
% The above example shows that even for moderate delays and minimal real
% parts accuracy is difficult to guarantee. The function
% |dde_stst_eig_cheb| gradually increases the degree of the polynomials
% until the discretization matrix reaches |maxsize| or until all dteceted
% eigenvalues with real part greater than |minimal_real_part| have residual
% less than |root_accuracy|. The functions |dde_stst_eig_bdf| and
% |dde_stst_eig_max| use heuristics (described by Verheyden et al) to
% determine required matrix size. All approaches will fail in some examples
% such that one may have to increase |minimal_real_part| closer to zero.
%% An extreme example
% Consider the linear DDE
% 
% $$x'(t)=-x(t)+y(t-\tau)\mbox{,\quad} y'(t)=-y(t)\mbox{.}$$
%
% This DDE has only one finite eigenvalue |-1| of double algebraic
% multiplicity. All other eigenvalues have disappeared at minus infinity,
% indpendent of the delay |tau|. For small delays |tau| all methods
% identify the spectrum correctly, but for larger delays spurious
% eigenvalues occupy part of the complex plane |-1 < Re z < 0| and small
% (in modulus) imaginary parts. This problem can not be remedied by
% refining the discretization.
A=cat(3,-eye(2),diag(1,1))
tausel=[20,50,40];
disc={'cheb','bdf','mxo'};
nmth=length(disc);
ntau=length(tausel);
lms_tstep=1e-3;
lw={'linewidth',2};
clear stability
for i=1:3
    mth=getfield(df_mthod('stst',disc{i}),'stability');
    % set method parameters such that matrix is large and
    % all eigenvalues of the discretized system are reported
    mth.minimal_real_part=-1.2;
    mth.minimal_time_step=lms_tstep;
    mth.maximal_time_step=lms_tstep;
    mth.max_number_of_eigenvalues=Inf;
    mth.max_newton_iterations=0;
    mth.remove_unconverged_roots=0;
    mth.maxsize=2050;
    mth.discard_accuracy_factor=Inf;
    stability{i}=feval(['dde_stst_eig_',disc{i}],A,tausel(i),...
        'method',mth);
    [~,ix]=sort(abs(stability{i}.l0+1),'ascend');
    figure(i);clf;hold on;
    pc2=plot(real(stability{i}.l0(ix(1:2))),imag(stability{i}.l0(ix(1:2))),...
        'bo',lw{:},'DisplayName','2 ev closest to -1');
    po=plot(real(stability{i}.l0(ix(3:end))),imag(stability{i}.l0(ix(3:end))),...
        'rx',lw{:},'DisplayName','spurious');
    title(sprintf('discretization %s, tau=%d',disc{i},tausel(i)));
    xlabel('Re');
    ylabel('Im');
    grid on
    legend([pc2,po]);
    hold off
    set(gca,lw{:},'fontsize',18,'fontname','courier','fontweight','bold');
end
%%
% We observe that the delay at which the problem becomes visible depends on
% the method. However, once we observe spurious eigenvalues near the
% correct eigenvalues, refining the discretization does not improve the
% results.