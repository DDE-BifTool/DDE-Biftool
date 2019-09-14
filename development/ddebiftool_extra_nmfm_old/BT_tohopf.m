function hopf=BT_tohopf(funcs,btpoint,freqs) %#ok<INUSD,INUSL>
hopf=btpoint;
hopf.kind='hopf';
hopf.omega=0;
hopf.v=btpoint.q0;
hopf=rmfield(hopf,'q0');
hopf=rmfield(hopf,'q1');
if isfield(hopf,'p0');
    hopf=rmfield(hopf,'p0');
end
if isfield(hopf,'p1');
    hopf=rmfield(hopf,'p1');
end
if isfield(btpoint,'stability')
    hopf.stability=btpoint.stability;
end
end 
  