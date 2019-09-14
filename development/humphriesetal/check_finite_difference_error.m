%% Compare normal forms obtained with finite diffferencing to those with symbolic derivatives
%
clear
s0=load('humphries_equilibria_symbolic=0.mat');
s1=load('humphries_equilibria_symbolic=1.mat');
%% Hopf Lyapunov coefficients
for i=1:length(s0.nf)
    fprintf('max of L1 error along Hopf curve %d is %g\n',...
        i,max(abs(s0.nf{i}-s1.nf{i})));
    fprintf('max of L1 error estimate along Hopf curve %d is %g\n\n',...
        i,max(abs(s0.nf{i}-s0.nflow{i})));
end
%% Hopf-Hopf interactions
for i=1:length(s0.nhh)
    names=fieldnames(s0.nhh(i).nmfm);
    err=struct2array(s0.nhh(i).nmfm)-struct2array(s1.nhh(i).nmfm);
    [abserr,ierr]=max(abs(err));
    errest=struct2array(s0.nhh(i).nmfm)-struct2array(s0.nhhlow(i).nmfm);
    [abserrest,ierrest]=max(abs(errest));
    fprintf('max of error in Hopf-Hopf %d is in %s: %g\n',...
        i,names{ierr},abserr);
    fprintf('max of error estimate in Hopf-Hopf %d is in %s: %g\n\n',...
        i,names{ierrest},abserrest);
end
    
    
    