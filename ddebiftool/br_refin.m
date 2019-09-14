function branch=br_refin(funcs,branch,x_m,y_m,axs)
%% add points to branch by clicking left mouse button,quit with middle or right mouse button
% function refined_branch=br_refin(funcs,branch,x_measure,y_measure,axs)
% INPUT:
%   funcs problem functions
%	branch branch of points
%	x_measure scalar measure for x-coordinate
%	y_measure scalar measure for y-coordinate
%	axs optional viewing axis
% OUTPUT:
%	refined_branch branch with refinements
% COMMENT:
%       add points to branch by clicking left mouse button
%       quit with middle or right mouse button

% (c) DDE-BIFTOOL v. 1.02, 21/09/2001
%
% $Id: br_refin.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
free_par=branch.parameter.free;

ll=length(branch.point);

if ll<=1
  err=ll;
  error('BR_REFIN: branch is incomplete, has only %d points.',err);
end;

stability=0;
replot=0;

if ~isempty(x_m)
  [ffx,llx]=br_measr(branch,x_m);
  if length(x_m.field)==9
    if strcmp(x_m.field,'stability')
      stability=1;
    end;
  end;
else
  ffx=(1:ll)';
  llx=ones(1,ll);
  replot=1;
end;
if ~isempty(y_m)
  [ffy,lly]=br_measr(branch,y_m);
  if length(y_m.field)==9 
    if strcmp(y_m.field,'stability')
      stability=1;
    end;
  end;
else
  ffy=(1:ll)';
  lly=ones(1,ll);
  replot=1;
end;

clf;
subplot(2,1,1);
br_plot(branch,x_m,y_m,'b');
hold on;
for i=1:ll
  if llx(i) && lly(i)
    plot(ffx(i,1:llx(i)),ffy(i,1:lly(i)),'b.');
  end;
end;
if exist('axs','var')
  axis(axs);
end;
ax1=axis;
subplot(2,1,2);
br_plot(branch,x_m,y_m,'b');
if exist('axs','var')
  axis(axs);
end;
ax2=axis;

[x,y,button]=ginput(1);

found=0;

while (button~=3 && button~=2)
  if button==1
    i=1;
    while i<ll
      if size(ffx,2)==1
        if (x<=ffx(i) && x>=ffx(i+1)) || (x<=ffx(i+1) && x>=ffx(i))
          if size(ffy,2)==1
            if (y<=ffy(i) && y>=ffy(i+1)) || (y<=ffy(i+1) && y>=ffy(i))
              found=1;
            end;
          else
            found=1;
          end;
        end;
      elseif size(ffy,2)==1
        if (y<=ffy(i) && y>=ffy(i+1)) || (y<=ffy(i+1) && y>=ffy(i))
          found=1;
        end;
      end;
      if found
        new_success=0;
        if stability 
          if isempty(branch.point(i).stability)
            branch.point(i).stability=p_stabil(funcs,branch.point(i),branch.method.stability);
            new_point=branch.point(i);
            new_success=1;
          elseif isempty(branch.point(i+1).stability)
            branch.point(i+1).stability=p_stabil(funcs,branch.point(i+1),branch.method.stability);
            new_point=branch.point(i+1); 
            new_success=1;
          end;
        end;
        if ~new_success
          new_point=p_axpy(1.0,branch.point(i),branch.point(i+1));
          new_point=p_axpy(0.5,new_point,[]);
          step_cnd=p_axpy(-1,branch.point(i),branch.point(i+1));
          [new_point,new_success]=p_correc(funcs,new_point,free_par, ...
						step_cnd,branch.method.point,i+1);
          if new_success
            if stability
              new_point.stability=p_stabil(funcs,new_point,branch.method.stability);
            else
              if isfield(branch.point(i+1),'stability');
                new_point.stability=[];
              end;
            end;
            branch.point(i+2:ll+1)=branch.point(i+1:ll);
            branch.point(i+1)=new_point;
          end;
        end;
        if new_success
          if ~isempty(x_m)
            fx=p_measur(new_point,x_m);
          else
            fx=i+1;
            ffx(i+1:ll)=ffx(i+1:ll)+1;
          end;
          if ~isempty(y_m)
            fy=p_measur(new_point,y_m);
          else
            fy=i+1;
            ffy(i+1:ll)=ffy(i+1:ll)+1;
          end;
          subplot(2,1,1);
          ffx(i+1:ll+1,:)=ffx(i:ll,:);
          llx(i+1:ll+1)=llx(i:ll);
          ffy(i+1:ll+1,:)=ffy(i:ll,:);
          lly(i+1:ll+1)=lly(i:ll);
          ffx(i+1,1:size(fx,2))=fx;
          llx(i+1)=length(fx);
          if size(fy,1)>1
            fy=fy';
          end;
          ffy(i+1,1:size(fy,2))=fy;
          lly(i+1)=length(fy);
          lx=min(llx(i),llx(i+1));
          ly=min(lly(i),lly(i+1));
          if replot
            hold off;
            br_plot(branch,x_m,y_m,'b');
            br_plot(branch,x_m,y_m,'b.');
          else
            if lx>1 && ly>1
              lx=min(lx,ly);
              ly=min(lx,ly);
            end; 
            plot(ffx(i:i+1,1:lx),ffy(i:i+1,1:ly),'r.')
            plot(ffx(i:i+1,1:lx),ffy(i:i+1,1:ly),'r-')
            lx=min(llx(i+1),llx(i+2));
            ly=min(lly(i+1),lly(i+2));
            plot(ffx(i+1:i+2,1:lx),ffy(i+1:i+2,1:ly),'r.')
            plot(ffx(i+1:i+2,1:lx),ffy(i+1:i+2,1:ly),'r-')
          end;
          ll=ll+1;
          i=i+1;
          subplot(2,1,2);
          hold off;
          br_plot(branch,x_m,y_m,'b');
        else
          disp('BR_REFIN warning: correction failed!');  
        end;
        found=0;
      end;
      i=i+1;
    end;
    subplot(2,1,1);
    axis(ax1);
    subplot(2,1,2);
    axis(ax2);
%  elseif button==2
%    if x<ax1(1) | y<ax1(3) | x>ax1(2) | y>ax1(4)
%      keyboard;
%    else
%      a=axis;
%      dx=(a(2)-a(1))/4;
%      dy=(a(4)-a(3))/4;
%      axis([x-dx x+dx y-dy y+dy]);
%    end;
  end;
  [x,y,button]=ginput(1);
  if (button==3 || button==2)
    subplot(2,1,1);
    ax=axis;
    if sum(ax~=ax1) 
      axis(ax1);
      button=0;
    end;
    subplot(2,1,2);
    ax=axis;
    if sum(ax~=ax2) 
      axis(ax2);
      button=0;
    end;
  end;
end;

subplot(2,1,1);
axis(ax1);
subplot(2,1,2);
axis(ax2);

return;


