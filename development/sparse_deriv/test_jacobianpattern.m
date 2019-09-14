%% test NodeColors
clear
Nrange=5:50;
sdim=3;
fail=false;
for k=1:length(Nrange)
    nbh1=ConnectionGraphFromGrid([Nrange(k),Nrange(k)],...
        'perbc',[true,false],'dist',[2,2]);
    nbh2=ConnectionsDouble(nbh1);
    clr=NodeColors(nbh2);
    Jp=JacobianPattern(clr,nbh1,'dim',sdim);
    fprintf('k=%d, N=%d\n',k,Nrange(k));
    for i=1:max(Jp.color)
        n1=Jp.rows(:,Jp.color==i);
        n2=sort(n1(:));
        ind=find(diff(n2)==0,1,'first');
        fprintf('clr=%d\n',i);
        if ~isempty(ind)
            fprintf('ind=%d fails\n',ind);
            fail=true;
            break
        end
    end
    if fail
        break;
    end
end
