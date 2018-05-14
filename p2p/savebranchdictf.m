%% 
% savebranchdictf(dir,fname,dictfname)
% save one of the folder branches to a python dictionary
%
function savebranchdictf(dir,pt,dictfname)
p=loadp(dir,pt);
savebranchdict(p,dictfname);
end