%% Simplified version of function_handle from OctSymPy
%
% returns only code string
%
% %Id$
%
% adapted from original functiona_handle (Copyright (C) 2014-2016 Colin
% B. Macdonald)
function code = dde_octave_code(f,fcnname,inputs)
  %% command as in function_handle, however with pre-set arguments
  cmd = { '(expr,fcnname,in_vars) = _ins' ...
            'from sympy.utilities.codegen import codegen' ...
            'try:' ...
           ['    out = codegen((fcnname,expr), "' 'octave' ...
            '", fcnname, header=1' ...
            ', argument_sequence=in_vars)'] ...
            'except ValueError as e:' ...
            '    return (False, str(e))' ...
            'return (True, out)' };
  
    [worked, out] = python_cmd (cmd, f, fcnname, inputs);

    if (~worked)
      if (strcmp(out, 'Language ''octave'' is not supported.'))
        error('function_handle: your SymPy has no octave codegen');
      else
        out
        error('function_handle: Some other error from SymPy code gen?  file a bug!');
      end
    end
    code = out{1}{2};

end


%!shared x,y,z
%! syms x y z

%!test
%! % basic test
%! h = function_handle(2*x);
%! assert(isa(h, 'function_handle'))
%! assert(h(3)==6)

%!test
%! % autodetect inputs
%! h = function_handle(2*x*y, x+y);
%! [t1, t2] = h(3,5);
%! assert(t1 == 30 && t2 == 8)

%!test
%! % specified inputs
%! h = function_handle(2*x*y, 'vars', [x y]);
%! assert(h(3,5)==30)
%! h = function_handle(2*x*y, x+y, 'vars', [x y]);
%! [t1, t2] = h(3,5);
%! assert(t1 == 30 && t2 == 8)

%!test
%! % cell arrays for vars list
%! h = function_handle(2*x*y, x+y, 'vars', {x y});
%! [t1, t2] = h(3,5);
%! assert(t1 == 30 && t2 == 8)
%! h = function_handle(2*x*y, x+y, 'vars', {'x' 'y'});
%! [t1, t2] = h(3,5);
%! assert(t1 == 30 && t2 == 8)

%!test
%! % cell arrays specfies order, overriding symvar order
%! h = function_handle(x*y, 12/y, 'vars', {y x});
%! [t1, t2] = h(3, 6);
%! assert(t1 == 18 && t2 == 4)
%! h = function_handle(x*y, 12/y, 'vars', [y x]);
%! [t1, t2] = h(3, 6);
%! assert(t1 == 18 && t2 == 4)

%!test
%! % cell arrays specfies order, overriding symvar order
%! h = function_handle(x*y, 12/y, 'vars', {y x});
%! [t1, t2] = h(3, 6);
%! assert(t1 == 18 && t2 == 4)
%! h = function_handle(x*y, 12/y, 'vars', [y x]);
%! [t1, t2] = h(3, 6);
%! assert(t1 == 18 && t2 == 4)

%!test
%! % Functions with different names in Sympy.
%! f = abs(x);  % becomes Abs(x)
%! h = function_handle(f);
%! assert(h(-10) == 10)
%! f = ceil(x);
%! h = function_handle(f);
%! assert(h(10.1) == 11)

%!test
%! % 'file' with empty filename returns handle
%! h = function_handle(2*x*y, 'file', '');
%! assert(isa(h, 'function_handle'))
%! assert(h(3,5)==30)
%! h = function_handle(2*x*y, 'vars', {x y}, 'file', '');
%! assert(isa(h, 'function_handle'))
%! assert(h(3,5)==30)

%!test
%! % output to disk
%! fprintf('\n')
%! if (exist ('OCTAVE_VERSION', 'builtin'))
%!   temp_file = tempname('', 'oct_');
%! else
%!   temp_file = tempname();
%! end
%! % allow loading function from temp_file
%! [temp_path, ans, ans] = fileparts(temp_file);
%! addpath(temp_path);
%! f = function_handle(2*x*y, 2^x, 'vars', {x y z}, 'file', temp_file);
%! assert( isa(f, 'function_handle'))
%! addpath(temp_path);  % Matlab 2014a needs this?
%! [a,b] = f(10,20,30);
%! assert (isnumeric (a) && isnumeric (b))
%! assert (a == 400)
%! assert (b == 1024)
%! if (exist ('OCTAVE_VERSION', 'builtin'))
%!   assert (unlink([temp_file '.m']) == 0)
%! else
%!   delete ([temp_file '.m'])
%! end
%! % remove temp_path from load path
%! rmpath(temp_path);

%!test
%! % output to disk: also works with .m specified
%! if (exist ('OCTAVE_VERSION', 'builtin'))
%!   temp_file = [tempname('', 'oct_') '.m'];
%! else
%!   temp_file = [tempname() '.m'];
%! end
%! % allow loading function from temp_file
%! [temp_path, ans, ans] = fileparts(temp_file);
%! addpath(temp_path);
%! f = function_handle(2*x*y, 2^x, 'vars', {x y z}, 'file', temp_file);
%! assert( isa(f, 'function_handle'))
%! addpath(temp_path);  % Matlab 2014a needs this?
%! [a,b] = f(10,20,30);
%! assert (isnumeric (a) && isnumeric (b))
%! assert (a == 400)
%! assert (b == 1024)
%! if (exist ('OCTAVE_VERSION', 'builtin'))
%!   assert (unlink(temp_file) == 0)
%! else
%!   delete (temp_file)
%! end
%! % remove temp_path from load path
%! rmpath(temp_path);

%!test
%! % non-scalar outputs
%! H = [x y z];
%! M = [x y; z 16];
%! V = [x;y;z];
%! h = function_handle(H, M, V);
%! [t1,t2,t3] = h(1,2,3);
%! assert(isequal(t1, [1 2 3]))
%! assert(isequal(t2, [1 2; 3 16]))
%! assert(isequal(t3, [1;2;3]))

%!test
%! % non-scalar outputs in .m files
%! H = [x y z];
%! M = [x y; z 16];
%! V = [x;y;z];
%! if (exist ('OCTAVE_VERSION', 'builtin'))
%!   temp_file = tempname('', 'oct_');
%! else
%!   temp_file = tempname();
%! end
%! % allow loading function from temp_file
%! [temp_path, ans, ans] = fileparts(temp_file);
%! addpath(temp_path);
%! h = function_handle(H, M, V, 'vars', {x y z}, 'file', temp_file);
%! assert( isa(h, 'function_handle'))
%! addpath(temp_path);  % Matlab 2014a needs this?
%! [t1,t2,t3] = h(1,2,3);
%! assert(isequal(t1, [1 2 3]))
%! assert(isequal(t2, [1 2; 3 16]))
%! assert(isequal(t3, [1;2;3]))
%! if (exist ('OCTAVE_VERSION', 'builtin'))
%!   assert (unlink([temp_file '.m']) == 0)
%! else
%!   delete ([temp_file '.m'])
%! end
%! % remove temp_path from load path
%! rmpath(temp_path);

%!test
%! % order of outputs is lexiographic
%! syms a A x y
%! f = y + 10*a + 100*x + 1000*A;
%! h = function_handle(f);
%! assert (h(1, 2, 3, 4) == 1000 + 20 + 300 + 4)
