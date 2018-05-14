%%
%  plotaf(dir,pt,fig): use plotsol_sp to plot the solution of first two 
%  components in 2D 
%  
function plotaf(dir,pt,fig)
figureHandle = figure(fig);clf;
subplot(1,2,1);
plotsol_sp(dir,pt,'pcmp',1,'sp');title('b_2');
d=load('vegcm.mat');
colormap(d.mycmap);
% Extract axes handles of all subplots from the figure
axesHandles = findobj(get(figureHandle,'Children'), 'flat','Type','axes');
% Set the axis property to square
axis(axesHandles,'image')
axis tight;
end
