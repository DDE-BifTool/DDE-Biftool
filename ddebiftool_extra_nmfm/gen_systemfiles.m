function s = gen_systemfiles(f,n,p,free_pars,taus,sys_dir,...
    vectorized,generate_multiplinearforms)
%% Generate the rhs file, the standard derivative file,
%  the vectorized standard derivative file,
%  the multilinear forms needed for calculating
%  the critical and parameter-dependent normal form coefficients
%  Input:
%    f the system
%    n number of components
%    p number of normal parameters
%    free_pars parameters used for unfolding codim 2 bifurcation
%    delay_ind array with indeces of nonzero delays
%    sys_dir directory where the system files are stored
%    vectorized boolean for generating vectorized standards derivative file
%    generate_multiplinearforms boolean for generating multilinear forms
%   Ouput:
%       s 1 success, 0 failed
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: gen_systemfiles.m 309 2018-10-28 19:02:42Z jansieber $
%
%% set defaults
if exist ('OCTAVE_VERSION', 'builtin') > 0
    gen_systemfiles_octave(f,n,p,free_pars,taus,sys_dir,...
        vectorized,generate_multiplinearforms);
else
m = length(taus)+1;
par = sym('par', [1,p+m-1]);
xx = sym('xx', [n,m]);
v = sym('v', [1,n]);

%% create directory
s=mkdir(sys_dir);
if ~s
    error('Could not create directory')
end
%% Generate rhs file
disp('Generating system file')

str='function f = sys_rhs(xx,par)\n\n';
for i=1:n
    str=strcat(str,sprintf('f(%d,:) = ',i),char(f(i)),';\n');
end
str=strcat(str,'\nend\n');

expression = 'par(\d+)';
replace = 'par($1)';
str = regexprep(str,expression,replace);

expression = 'xx(\d+)_(\d+)';
replace = 'xx($1,$2,:)';
str = regexprep(str,expression,replace);

expression = '*';
replace = '.*';
str = regexprep(str,expression,replace);

expression = '\^';
replace = '.^';
str = regexprep(str,expression,replace);

% save to file
fid = fopen(strcat(sys_dir,'sys_rhs.m'), 'w');
fprintf(fid, str);
fclose(fid);
%% Generate standard derivative file
disp('Generating standard derivative file')

str = strcat(['function J = sys_deri(xx,par,nx,np,v)\n\n'...
    'J = [];\n\n'...
    'if length(nx) == 1 && isempty(np) && isempty(v)\n'...
    '\tswitch nx\n']);
for nx=0:m-1
    str = strcat(str,sprintf('\tcase %d',nx),...
        '\n\t\tJ = [');
    for i=1:n
        for j=1:n
            deri = diff(f(i), xx(j,nx+1));
            str = strcat(str,char(deri));
            if j < n
                str = strcat(str,', ');
            end
        end
        if i < n
            str = strcat(str,'; ');
        end
    end
    str = strcat(str,'];\n');
end
str = strcat(str,['\tend\n'...
    'elseif isempty(nx) && length(np) == 1 && isempty(v)\n'...
    '\tswitch np\n']);
for np=1:p
    str = strcat(str,sprintf('\tcase %d',np),...
        '\n\t\tJ = [');
    for i=1:n
        deri = diff(f(i), par(np));
        str = strcat(str,char(deri));
        if i < n
            str = strcat(str,'; ');
        end
    end
    str = strcat(str,'];\n');
end
str = strcat(str,['\tend\n'...
    'elseif length(nx) == 1 && length(np) == 1 && isempty(v)\n'...
    '\tswitch nx\n']);
for nx=0:m-1
    str = strcat(str,sprintf('\tcase %d',nx),...
        '\n\t\tswitch np\n');
    for np=1:p
        str = strcat(str,sprintf('\t\tcase %d',np),...
            '\n\t\t\tJ = [');
        for i=1:n
            for j=1:n
                deri = diff(f(i), par(np), xx(j,nx+1));
                str = strcat(str,char(deri));
                if j < n
                    str = strcat(str,', ');
                end
            end
            if i < n
                str = strcat(str,'; ');
            end
        end
        str = strcat(str,'];\n');
    end
    str = strcat(str,'\t\tend\n');
end
str = strcat(str,['\tend\n'...
    'elseif length(nx) == 2 && isempty(np) && ~isempty(v)\n'...
    '\tnx1 = nx(1); nx2 = nx(2);\n'...
    '\tswitch nx1\n']);
for nx1=0:m-1
    str = strcat(str,sprintf('\tcase %d',nx1),...
        '\n\t\tswitch nx2\n');
    for nx2=0:m-1
        str = strcat(str,sprintf('\t\tcase %d',nx2),...
            '\n\t\t\tJ = [');
        for i=1:n
            for j=1:n
                deri1=0;
                for s=1:n
                    deri1=deri1+diff(f(i),xx(s,nx1+1))*v(s);
                end
                deri = diff(deri1, xx(j,nx2+1));
                str = strcat(str,char(deri));
                if j < n
                    str = strcat(str,', ');
                end
            end
            if i < n
                str = strcat(str,'; ');
            end
        end
        str = strcat(str,'];\n');
    end
    str = strcat(str,'\t\tend\n');
end
str = strcat(str,['\tend\n'...
    'end\n'...
    'if isempty(J)\n\tdisplay([nx np size(v)]);\n'...
    '\terror(''SYS_DERI: requested derivative could not be computed!'');\nend']);

expression = 'par(\d+)';
replace = 'par($1)';
str = regexprep(str,expression,replace);

expression = 'xx(\d+)_(\d+)';
replace = 'xx($1,$2)';
str = regexprep(str,expression,replace);

expression = 'v(\d+)';
replace = 'v($1)';
str = regexprep(str,expression,replace);

% save to file
fid = fopen(strcat(sys_dir,'sys_deri.m'), 'w');
fprintf(fid, str);
fclose(fid);
%% Generate vectorized standard derivative file
if vectorized
    disp('Generating vectorized standard derivative file')
    str = strcat(['function J = sys_deri_vec(xx,par,nx,np,v)\n\n'...
        'J = [];\n\n'...
        'I=ones(size(xx(1,1,:)));\n\n'...
        'if length(nx) == 1 && isempty(np) && isempty(v)\n'...
        '\tswitch nx\n']);
    for nx=0:m-1
        str = strcat(str,sprintf('\tcase %d',nx),...
            '\n\t\tJ = [');
        for i=1:n
            for j=1:n
                deri = diff(f(i), xx(j,nx+1));
                str = strcat(str,'(',char(deri),'*I)');
                if j < n
                    str = strcat(str,', ');
                end
            end
            if i < n
                str = strcat(str,'; ');
            end
        end
        str = strcat(str,'];\n');
    end
    str = strcat(str,['\tend\n'...
        'elseif isempty(nx) && length(np) == 1 && isempty(v)\n'...
        '\tswitch np\n']);
    for np=1:p
        str = strcat(str,sprintf('\tcase %d',np),...
            '\n\t\tJ = [');
        for i=1:n
            deri = diff(f(i), par(np));
            str = strcat(str,'(',char(deri),'*I)');
            if i < n
                str = strcat(str,'; ');
            end
        end
        str = strcat(str,'];\n');
    end
    str = strcat(str,['\tend\n'...
        'elseif length(nx) == 1 && length(np) == 1 && isempty(v)\n'...
        '\tswitch nx\n']);
    for nx=0:m-1
        str = strcat(str,sprintf('\tcase %d',nx),...
            '\n\t\tswitch np\n');
        for np=1:p
            str = strcat(str,sprintf('\t\tcase %d',np),...
                '\n\t\t\tJ = [');
            for i=1:n
                for j=1:n
                    deri = diff(f(i), par(np), xx(j,nx+1));
                    str = strcat(str,'(',char(deri),'*I)');
                    if j < n
                        str = strcat(str,', ');
                    end
                end
                if i < n
                    str = strcat(str,'; ');
                end
            end
            str = strcat(str,'];\n');
        end
        str = strcat(str,'\t\tend\n');
    end
    str = strcat(str,['\tend\n'...
        'elseif length(nx) == 2 && isempty(np) && ~isempty(v)\n'...
        '\tnx1 = nx(1); nx2 = nx(2);\n'...
        '\tswitch nx1\n']);
    for nx1=0:m-1
        str = strcat(str,sprintf('\tcase %d',nx1),...
            '\n\t\tswitch nx2\n');
        for nx2=0:m-1
            str = strcat(str,sprintf('\t\tcase %d',nx2),...
                '\n\t\t\tJ = [');
            for i=1:n
                for j=1:n
                    deri1=0;
                    for s=1:n
                        deri1=deri1+diff(f(i),xx(s,nx1+1))*v(s);
                    end
                    deri = diff(deri1, xx(j,nx2+1));
                    str = strcat(str,'(',char(deri),'*I)');
                    if j < n
                        str = strcat(str,', ');
                    end
                end
                if i < n
                    str = strcat(str,'; ');
                end
            end
            str = strcat(str,'];\n');
        end
        str = strcat(str,'\t\tend\n');
    end
    str = strcat(str,['\tend\n'...
        'end\n'...
        'if isempty(J)\n\tdisplay([nx np size(v)]);\n'...
        '\terror(''SYS_DERI: requested derivative could not be computed!'');\nend']);
    
    expression = 'par(\d+)';
    replace = 'par($1)';
    str = regexprep(str,expression,replace);
    
    expression = 'xx(\d+)_(\d+)';
    replace = 'xx($1,$2,:)';
    str = regexprep(str,expression,replace);
    
    expression = 'v(\d+)';
    replace = 'v($1,:,:)';
    str = regexprep(str,expression,replace);
    
    expression = '*';
    replace = '.*';
    str = regexprep(str,expression,replace);
    
    expression = '\^';
    replace = '.^';
    str = regexprep(str,expression,replace);
    
    % save to file
    fid = fopen(strcat(sys_dir,'sys_deri_vec.m'), 'w');
    fprintf(fid, str);
    fclose(fid);
end

%% Generate free parameter file
if generate_multiplinearforms
    disp('Generating parameter file file')
    fid = fopen(strcat(sys_dir,'get_free_pars.m'), 'w');
    fprintf(fid, 'function free_pars = get_free_pars();\r\n');
    fprintf(fid, 'free_pars=%s;\r\n',mat2str(free_pars));
    fprintf(fid, 'end');
    fclose(fid);
end
%% Generate sys_tau file
disp('Generating sys_tau file')
fid = fopen(strcat(sys_dir,'sys_tau.m'), 'w');
fprintf(fid, 'function tau = sys_tau();\r\n');
fprintf(fid, 'tau=%s;\r\n',mat2str(taus));
fprintf(fid, 'end');
fclose(fid);

%% Generate multilinear form files
if generate_multiplinearforms
    disp('Generating multilinear form files for normal form computation')
    %% B
    phi0 = sym('phi0', [n*m,1]);
    phi1 = sym('phi1', [n*m,1]);
    
    Bphi1phi1=jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1;
    
    matlabFunction(Bphi1phi1, 'file',strcat(sys_dir,'sys_B'),...
        'vars',{xx,par,phi0,phi1});
    
    %% B2
    phi0 = sym('phi0', [n*m,1]);
    phi1 = sym('phi1', [n*m,1]);
    psi0=sym('psi0',[2,1]);
    
    Bphi1phi1psi0=jacobian(jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1,...
        [par(free_pars(1)); par(free_pars(2))])*psi0;
    
    matlabFunction(Bphi1phi1psi0, 'file',strcat(sys_dir,'sys_B2'),...
        'vars',{xx,par,phi0,phi1,psi0});
    
    %% C
    phi0 = sym('phi0', [n*m,1]);
    phi1 = sym('phi1', [n*m,1]);
    phi2 = sym('phi2', [n*m,1]);
    Bphi1phi1phi2=jacobian(jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1,xx(:))*phi2;
    matlabFunction(Bphi1phi1phi2, 'file',strcat(sys_dir,'sys_C'),...
        'vars',{xx,par,phi0,phi1,phi2});
    
    %% J1
    J1=jacobian(f,[par(free_pars(1)); par(free_pars(2))]);
    matlabFunction(J1, 'file',strcat(sys_dir,'sys_J1'), 'vars',{xx,par});
    
    %% J2
    psi0=sym('psi0',[2,1]);
    psi1=sym('psi1',[2,1]);
    J2=jacobian(jacobian(f,[par(free_pars(1)); par(free_pars(2))])*psi0,...
        [par(free_pars(1)); par(free_pars(2))])*psi1;
    matlabFunction(J2, 'file',strcat(sys_dir,'sys_J2'), 'vars',{xx,par,psi0,psi1});
    
    %% A1
    phi0=sym('phi0', [n*m,1]);
    psi1=sym('psi0',[2,1]);
    A1=jacobian(jacobian(f,xx(:))*phi0,[par(free_pars(1)); par(free_pars(2))])*psi1;
    matlabFunction(A1, 'file',strcat(sys_dir,'sys_A1'),...
        'vars',{xx,par,phi0,psi0});
    
    %% D
    phi0 = sym('phi0', [n*m,1]);
    phi1 = sym('phi1', [n*m,1]);
    phi2 = sym('phi2', [n*m,1]);
    phi3 = sym('phi3', [n*m,1]);
    Dphi1phi1phi2phi3=jacobian(jacobian(jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1,xx(:))*phi2,xx(:))*phi3;
    matlabFunction(Dphi1phi1phi2phi3, 'file',strcat(sys_dir,'sys_D'),...
        'vars',{xx,par,phi0,phi1,phi2,phi3});
    
    %% E
    phi0 = sym('phi0', [n*m,1]);
    phi1 = sym('phi1', [n*m,1]);
    phi2 = sym('phi2', [n*m,1]);
    phi3 = sym('phi3', [n*m,1]);
    phi4 = sym('phi4', [n*m,1]);
    Dphi1phi1phi2phi3phi4=jacobian(jacobian(jacobian(jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1,xx(:))*phi2,xx(:))*phi3,xx(:))*phi4;
    matlabFunction(Dphi1phi1phi2phi3phi4, 'file',strcat(sys_dir,'sys_E'),...
        'vars',{xx,par,phi0,phi1,phi2,phi3,phi4});
    
    %% C1
    phi0 = sym('phi0', [n*m,1]);
    phi1 = sym('phi1', [n*m,1]);
    phi2 = sym('phi2', [n*m,1]);
    psi1=sym('psi0',[2,1]);
    C1phi1phi1phi2psi1=jacobian(jacobian(jacobian(jacobian(f,xx(:))*phi0,xx(:))*phi1,xx(:))*phi2,[par(free_pars(1)); par(free_pars(2))])*psi1;
    matlabFunction(C1phi1phi1phi2psi1, 'file',strcat(sys_dir,'sys_C1'),...
        'vars',{xx,par,phi0,phi1,phi2,psi1});
    
    %% generate multilinear form file as used by the critical normal form computation
    str=['function y = sys_mfderi(xx,par,varargin)\n\n'...
        'if nargin == 2\n'...
        '\terror(''SYS_MFDERI: no arguments.'');\n'...
        'elseif nargin > 7\n'...
        '\terror(''SYS_MFDERI: too many arguments.'');\n'...
        'end\n\n'...
        'y=0;\n\n'...
        'numarg = nargin - 2;\n\n'];
    
    str = strcat(str,['switch numarg\n'...
        '\tcase 2\n'...
        '\t\tu1 = varargin{1}; u2 = varargin{2};\n'...
        '\t\ty=sys_B(xx,par,u1(:),u2(:));\n']);
    
    str = strcat(str, ['\tcase 3\n'...
        '\t\tu1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3};\n'...
        '\t\ty=sys_C(xx,par,u1(:),u2(:),u3(:));\n']);
    
    str = strcat(str, ['\tcase 4\n'...
        '\t\tu1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4};\n'...
        '\t\ty=sys_D(xx,par,u1(:),u2(:),u3(:),u4(:));\n']);
    
    str = strcat(str, ['\tcase 5\n'...
        '\t\tu1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4}; u5 = varargin{5};\n'...
        '\t\ty=sys_E(xx,par,u1(:),u2(:),u3(:),u4(:),u5(:));\n']);
    
    str=strcat(str,'\nend');
    
    % save to file
    fid = fopen(strcat(sys_dir,'sys_mfderi.m'), 'w');
    fprintf(fid, str);
    fclose(fid);
    %% done
end
disp('done')
end
s=1;
end
