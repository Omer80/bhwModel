%%
% plotaf2: plot 1st and 2nd component of p.u in struct p.
%  plotsol(p,varargin) 
%  
function plotaf2(p,fig)
figure(fig);clf;
d=load('vegcm.mat');
subplot(1,2,1);
plotsol_sp(p,'pcmp',1,'sp');title('b_2');
colormap(d.mycmap);
axis tight;
end
