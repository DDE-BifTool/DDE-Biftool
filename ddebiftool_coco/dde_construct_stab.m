function [prob,stabdata]=dde_construct_stab(prob,varargin)
%% constructor for dde-biftool stability
% $Id: dde_construct_stab.m 362 2019-07-14 15:49:40Z jansieber $
%%
default={'id','','bd_stability',true,'pointtypes',@pointtype_list};
[options,pass_on]=dde_set_options(default,varargin,'pass_on','data');
if isempty(options.id)
    fid='dde';
else
    fid=options.id;
end
id=options.id;
[ptdata,uidx] = coco_get_func_data(prob, fid,'data','uidx');
%% add stability monitors
tstfid=coco_get_id(id,'test');
pids = coco_get_id(id,{'USTAB','USTAB_C','USTAB_P','USTAB_N','CRIT'});
prob = coco_add_chart_data(prob, tstfid, [], []);
pt=dde_stabilitykind(ptdata.info.point,ptdata,false);
fcn=options.pointtypes('kind',pt.kind);
hopftrlabel={'HB','TR'};
isfloquet=fcn.stab(-2)>0;
hopftrlabel=hopftrlabel{isfloquet+1};
stabdata=struct('ptid',fid,'ownfid',tstfid,'id',id,'pids',{pids},'evfcn',fcn,...
    'args',{pass_on},'bd_stability',options.bd_stability,...
    'hopftrlabel',hopftrlabel);
prob=coco_add_func(prob,tstfid,@dde_stab,stabdata,'regular',...
    pids,'uidx',uidx, 'passChart');
%% set up event handling for bifurcations
prob=coco_add_event(prob, @evhan_stability, stabdata, 'SP', pids{2}, 0);
prob=coco_add_event(prob, @evhan_stability, stabdata, 'SP', pids{3}, 0);
prob=coco_add_event(prob, @evhan_stability, stabdata, 'SP', pids{4}, 0);
if options.bd_stability
    tstoutput=coco_get_id(tstfid,'out');
    prob = coco_add_slot(prob, tstoutput, @bddat, stabdata, 'bddat');
end
end
%%
function [data,chart, y]=dde_stab(prob,data,chart,u)
ptdata=coco_get_func_data(prob,data.ptid,'data');
cdata=coco_get_chart_data(chart,data.ownfid);
if ~isfield(cdata,'stability')
    spoint=dde_stabilitykind(u,ptdata,data.args);
    ev.stability=spoint.stability;
    ev.lambda=data.evfcn.getev(spoint);
    ev.trivial_ev=data.evfcn.triv(spoint);
    ev.stab=data.evfcn.stab;
    ev.isFloquet=ev.stab(-2)>0;
    chart=coco_set_chart_data(chart,data.ownfid,ev);
else
    ev=cdata;
end
y=lambda2tst(ev.lambda,ev.trivial_ev,ev.stab);
end
%%
function [data, cseg, msg] = evhan_stability(~, data, cseg, cmd, msg)
%% Bifurcation event handler fr bifurcations based on eigenvalues.
fid = data.ownfid;
switch cmd
  case 'init'
    if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
      msg.action = 'finish';
    elseif strcmp(msg.action, 'locate')
      msg.action = 'warn';
    else
      cdata = coco_get_chart_data(cseg.ptlist{1}, fid);
      tst0=lambda2tst(cdata.lambda,cdata.trivial_ev,cdata.stab);
      cdata = coco_get_chart_data(cseg.ptlist{end}, fid);
      tst1=lambda2tst(cdata.lambda,cdata.trivial_ev,cdata.stab);
      if tst0(2)~=tst1(2) && tst0(1)~=tst1(1) && mod(tst0(1)-tst1(1),2)==0 
          msg.point_type = data.hopftrlabel;
          msg.action     = 'locate';
      elseif mod(tst0(1)-tst1(1),2)==1 && tst0(3)~=tst1(3)
          msg.point_type = 'SN';
          msg.action     = 'locate';
      elseif mod(tst0(1)-tst1(1),2)==1 && tst0(4)~=tst1(4)
          msg.point_type = 'PD';
          msg.action     = 'locate';
      else
          msg.point_type = '';
          msg.action     = 'warn';
          msg.wmsg       = 'collision of positive eigenvalues';
      end
      msg.idx = 1;
    end
      case 'check'
      msg.action='add';
      msg.finish = true;
end
end
%%
function [y,non_triv]=lambda2tst(lambda_inp,trivial_ev,stab)
[lambda,~,non_triv]=dde_separate_complex(lambda_inp,trivial_ev);
y(1)=sum(stab(lambda)>=0);
y(2)=sum(stab(lambda)>=0&imag(lambda)>0);
y(3)=sum(stab(lambda)>=0&imag(lambda)==0&real(lambda)>=0);
y(4)=sum(stab(lambda)>=0&imag(lambda)==0&real(lambda)<=0);
if ~isempty(lambda)
    [~,ix]=min(abs(stab(lambda)));
    y(5)=real(lambda(ix))+1i*abs(imag(lambda(ix)));
else
    y(5)=NaN;
end
y(2:4)=mod(2*y(2:4)+1,4)-2;
end
%%
function [data, res] = bddat(prob, data, command, varargin) %#ok<INUSL>
% Append eigenvalues of Jacobian to BD.

res = {};
switch command
  case 'init'
    res = coco_get_id(data.id,'stability');
  case 'data'
    chart = varargin{1}; % Current chart
    cdata = coco_get_chart_data(chart, data.ownfid);
    if ~isempty(cdata) && isfield(cdata, 'stability')
      res = cdata.stability;
    end
end
end
