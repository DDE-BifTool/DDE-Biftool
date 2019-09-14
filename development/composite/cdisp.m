function cdisp(st,levels,indent,nind)
% cdisp Display a structure and its contents
% 20/1/01 Matlab5 W.Whiten with Bob Hsu & Andrew Wang
%
% cdisp(st)
% cdisp(st,levels,indent,nind)
% st A structure/array to be displayed
% levels Number of levels to display (optional)
% indent Characters to indent for this level (optional)
% nind Number of spaces to indent each level (optional)

% set default values
if(nargin<2)
   levels=5;
  end
if(nargin<3)
   indent='-';
  end
if(nargin<4)
   nind=3;
  end
mark=''; % preceeds type names
structindent='.'; % indent for struct names

% check if last level into structure
if(levels<=0)
   disp(st)
   return
  end

% size of input
nst=size(st);

% process the different types of input
switch(class(st))

  case('char')
   disp([indent, mark, 'char'])
   if(isempty(st))
      disp(' ''''')
     else
      disp(st)
     end

  case('double') % and logical
   if(islogical(st))
      if(isempty(st))
         if(length(nst)==2 & sum(nst)==0)
            disp([indent, mark, 'logical'])
            disp(' []')
           else
            disp([indent, mark, 'logical ', num2str(nst)])
            disp(' []')
           end
        else
         disp([indent, mark, 'logical'])
         disp(st)
      end

     else % must be double
      if(isempty(st))
         if(length(nst)==2 & sum(nst)==0)
            disp([indent, mark, 'double'])
            disp(' []')
           else
            disp([indent, mark, 'double ', num2str(nst)])
            disp(' []')
           end
        else
         disp([indent, mark, 'double'])
         disp(st)
      end
     end

  case('cell')
   if(isempty(st))
      if(length(nst)==2 & sum(nst)==0)
         disp([indent, mark, 'cell'])
         disp(' {}')
        else
         disp([indent, mark, 'cell ', num2str(nst)])
         disp(' {}')
        end

     elseif(prod(nst)==1) % one element
      disp([indent, mark, 'cell '])
      cdisp(st{1},levels-1,[indent, blanks(nind)],nind)
      
     else % multiple elements
      ist=ones(size(nst));
      i=1;
      while(i<=length(ist))
         disp([indent, mark, 'cell ', num2str(ist)])
         cdisp(st{subs1d(nst,ist)},levels-1, ...
            [indent, blanks(nind)],nind)

         i=1;
         while(i<=length(ist))
            ist(i)=ist(i)+1;
            if(ist(i)>nst(i))
               ist(i)=1;
              else
               break
              end

            i=i+1;

           end
        end
     end

  case('struct')
   if(isempty(st))
      if(length(nst)==2 & sum(nst)==0)
         disp([indent, mark, 'struct'])
         disp(fieldnames(st)')
         disp(' []')
        else
         disp([indent, mark, 'struct ', num2str(nst)])
         disp(fieldnames(st)')
         disp(' []')
        end

     elseif(prod(nst)==1) % one element
      fn=fieldnames(st);
      disp([indent, mark, 'struct '])
      indent1=[indent, structindent];
      for i=1:length(fn);
         disp([indent1, fn{i}])
         cdisp(getfield(st,fn{i}),levels-1, ...
            [indent, blanks(nind)],nind)
        end

     else % multiple elements

      ist=ones(size(nst));
      i=1;
      fn=fieldnames(st);
      indent1=[indent, structindent];
      while(i<=length(ist))
         disp([indent, mark, 'struct ', num2str(ist)])
         for i=1:length(fn);
            disp([indent1, fn{i}])
            cdisp(getfield(st,{subs1d(nst,ist)},fn{i}), ...
               levels-1,[indent, blanks(nind)],nind)
           end

         i=1;
         while(i<=length(ist))
            ist(i)=ist(i)+1;
            if(ist(i)>nst(i))
               ist(i)=1;
              else
               break
              end

            i=i+1;

           end
        end
     end

  otherwise
   disp([indent, mark, class(st)])
   ms=methods(class(st));
   msx=strmatch(ms,'cdisp','exact');
   if(isempty(msx))
      disp(st)
     else
      cdisp(st,levels-1,[indent1, blanks(nind)],nind)
     end

  end

return



function i=subs1d(dim,idx)
% subs1d Convert multi dimendsional subscript to linear (1d) subscript
% 3/2/01 Matlab5 W.Whiten
%
% i=subs1d(dim,idx)
% dim Dimension information for an array (from size)
% idx Row vector of subscripts for array
%
% i Value of a one dimensional subscript for use in array

dim=cumprod(dim);
i=sum([1,dim(1:end-1)].*(idx-1))+1;

return