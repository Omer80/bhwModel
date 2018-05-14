function savebranchdict(p,dictfname)
data=[];
data.ind=squeeze(p.branch(1,:)');
data.prec=squeeze(p.branch(4,:)');
data.Y2=squeeze(p.branch(end,:)');
data.Btot = squeeze(p.branch(end-4,:)');
save(dictfname,'-struct','data');
end