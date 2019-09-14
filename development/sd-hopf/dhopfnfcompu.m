% test ddebiftool normal form comp for double Hopf
% onbly works for constant delay DDE so we implelement our expanded
% truncated nonlinear constant delay DDE

clearvars; close all; tic;
diaryname=strcat('diary',date,'.txt');
diary(diaryname)

%addpath('C:/Program Files/MATLAB/dde_biftool_v3.1/ddebiftool/',...
%    'C:/Program Files/MATLAB/dde_biftool_v3.1/ddebiftool_extra_psol/',...
%    'C:/Program Files/MATLAB/dde_biftool_v3.1/ddebiftool_extra_nmfm/',...
%    'C:/Program Files/MATLAB/dde_biftool_v3.1/ddebiftool_utilities/');

funcs=set_funcs('sys_rhs', @(xx,par) sys_cub_rhs(xx,par),'sys_tau', sys_tau3 );  
    % use finite differences for derivatives

%% set parameters
% epsilon=par(1); gamma=par(2); kappa=par(3:4); a=par(5:6); c=par(7:8);
% a2=par(9:11);
% a2=(2a1 a1+a2 2a2);
% a3=par(12:15);
% a3=(3a1 2a1+a2 a1+2a2 3a2);
epsilon=1; gamma=4.75; a=[1.3 6]; c=[1 1]; kappa(2)=3; kappa(1)=gamma-kappa(2);
a2=[2*a(1) a(1)+a(2)  2*a(2)];
a3=[3*a(1) 2*a(1)+a(2) a(1)+2*a(2)  3*a(2)];
par0=[epsilon gamma kappa a c a2 a3];

% parameters a,a2,a3 define delays:
sys_tau3=@()[5:6 9:11 12:15];

%% continue equilibrium
% continue in kappa_1
bifpar=3; k1upper=14; % kappa_1 \in [0,14]

% set up branch
[eqbr,suc]=gen_stst(funcs,'contpar',bifpar,'x',0,'parameter',par0,...
    'max_bound',[bifpar,k1upper],'max_step',[bifpar,0.1]);
if ~suc
    error('equilibrium not found');
end

% continue branch
figure(1); eqbr=br_contn(funcs,eqbr,200);

disp('Computing steady state stability')
% find bifurcations
eqbr.method.stability.maximal_time_step=0.05;
eqbr.method.stability.minimal_time_step=0.0005;
eqbr.method.stability.minimal_real_part=-0.05;

%[eqnunst,dom,triv_defect,eqbr.point]=...
[~,~,~,eqbr.point]=...
    GetStability(eqbr,'funcs',funcs,'points',2:length(eqbr.point));

disp('Detecting steady state bifurcations')
[eqbrbifs,suc] = br_bifdet(funcs,eqbr);
if ~suc
    error('bifs not found');
end
fpi=br_getflags(eqbrbifs);
if numel(fpi) ~= numel(fpi(1,:))
  disp('warning bifs other than hopf found')
end

disp('Computing hopf curve')
%% we know from previous comps that for these parameters second branch is one with all the double Hopfs
start_ind=fpi(bif2num('hopf'),2);

[hopfbr1, suc] = SetupHopf(funcs,eqbrbifs, start_ind, 'contpar', [3 4], 'dir',4); %, 'step', cstepsize);
hopfbr1.parameter.min_bound=[3 0; 4 0];
hopfbr1.parameter.max_bound=[3 k1upper; 4 gamma];
hopfbr1.parameter.max_step=[3 0.1; 4 0.1];
hopfbr1.method.stability.minimal_time_step = 0.001; % default 0.01
%hopfbr1.method.stability.minimal_real_part=-0.2; %-0.5;
hopfbr1.method.bifurcation.minimal_real_part = -0.03;
hopfbr1.method.bifurcation.secant_tolerance=1e-5;
hopfbr1.method.bifurcation.correction_tolerance=1e-5;

figure(2); hopfbr1=br_contn(funcs,hopfbr1,200);    
hopfbr1=br_rvers(hopfbr1);
figure(2); hopfbr1=br_contn(funcs,hopfbr1,200);    

disp('Computing Hopf curve stability')
[eqnunst,dom,triv_defect,hopfbr1.point]=GetStability(hopfbr1,'funcs',funcs,'exclude_trivial',true);

figure(3);
k1=cell2mat(arrayfun(@(x)x.parameter(3)',hopfbr1.point,'uniformoutput',false));
omega=cell2mat(arrayfun(@(x)x.omega',hopfbr1.point,'uniformoutput',false));
plot(k1,omega)

disp('Detecting bifurcations from Hopf curve')

WherePerBifs=find(eqnunst(2:end)-eqnunst(1:end-1)~=0);
disp(['Found ',num2str(numel(WherePerBifs)),' bifurcations'])

for i=1:numel(WherePerBifs)
    disp(' ')
    disp(['Potential Bifurcation #',num2str(i)])    
    disp(['Between branch points ',num2str(WherePerBifs(i)),' and ',num2str(WherePerBifs(i)+1)])
    k12f=cell2mat(arrayfun(@(x)x.parameter(hopfbr1.parameter.free)',hopfbr1.point(WherePerBifs(i)),'uniformoutput',false));
    k12l=cell2mat(arrayfun(@(x)x.parameter(hopfbr1.parameter.free)',hopfbr1.point(WherePerBifs(i)+1),'uniformoutput',false));
    disp(['For ',num2str(k12f(1)),'<k1<',num2str(k12l(1)),' : ',num2str(k12f(2)),'<k2<',num2str(k12l(2))])
    disp(['Number of unstable eigenvalues changes from ',num2str(eqnunst(WherePerBifs(i))),' to ',num2str(eqnunst(WherePerBifs(i)+1))])
    
    if abs(eqnunst(WherePerBifs(i))-eqnunst(WherePerBifs(i)+1))==2
        disp('Refining potential double Hopf point')
        ttoc=toc;
        [nf,nflow,br_ref,indbif]=HopfHopfNormalform(funcs,hopfbr1,[WherePerBifs(i) WherePerBifs(i)+1]);
        tttoc=toc;
        disp(['This normal form computed in ',num2str(tttoc-ttoc),' seconds'])
        disp(' ')
        disp(['Double Hopf at k1=',num2str(nf.parameter(3),8),', k2=',num2str(nf.parameter(4),8)]);
        disp(['omega1=',num2str(nf.omega1,8),', omega2=',num2str(nf.omega2,8)]);
        disp('High Order Approx:')
        disp(nf.nmfm)
        disp('Low Order Approx:')
        disp(nflow.nmfm)
    end

end

%% br_bifdet fails to correct all the bifurcations
% [hopfbr1,suc] = br_bifdet(funcs,hopfbr1);


%% tidy up
%rmpath('C:/Program Files/MATLAB/dde_biftool_v3.1/ddebiftool/',...
%    'C:/Program Files/MATLAB/dde_biftool_v3.1/ddebiftool_extra_psol/',...
%    'C:/Program Files/MATLAB/dde_biftool_v3.1/ddebiftool_extra_psol/',...
%    'C:/Program Files/MATLAB/dde_biftool_v3.1/ddebiftool_utilities/');

tttoc=toc;
disp(['Total Runtime:',num2str(tttoc),' seconds'])
diary('off')


