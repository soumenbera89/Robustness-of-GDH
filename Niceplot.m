%% To create nice looking plot for publishing paper......
function niceplot(papersize, margin, frontsize)
% Function which produces nice looking plot..
% and set up the page for nice printing..
if nargin==0
    papersize = 12;
    margin = 0.5;
    frontsize = 8;
elseif nargin==1
   margin = 0.5;
   frontsize = 12;
else nargin==2
    frontsize = 12;
end
set(get(gca,'xlabel'),'FontSize', frontsize, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', frontsize, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', frontsize, 'FontWeight', 'Bold');
box off; axis square;
set(gca,'LineWidth',1.5);
set(gca,'FontSize',12);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [papersize papersize]);
set(gcf,'PaperPosition',[margin margin papersize-2*margin papersize-2*margin]);
set(gcf,'PaperPositionMode','Manual');
%set(gca,'XTick',0:4000);
%set(gca,'YTickLabelMode','manual')
%set(gca, 'YTick', 0:250);
%set(gca,'XTick',0:1);
%set(gca,'YTick',0:250);
%set(gca, 'XTickLabel', {'0';'10';'20';'30';'40';'50';'60';'70';'80'});
set(gca, 'XTickLabel', {'0';'2'});
%ax=gca;
%ax.XGrid = 'on';
%ax.YGrid = 'on';
%ax.LineWidth = 0.25;
%ax.XTickLabel={'0';'';'10';'';'20';'';'30';'';'40';'';'50';};
xlabel('Time(min)','FontWeight','Bold','FontSize',8);
ylabel('Glutamate(mM)','FontWeight','Bold','FontSize',8);
%ylabel('[GDH.Glu]^{-1}(mM)^{-1}','FontWeight','Bold','FontSize',12);
%xlabel('[KGA]^{-1}(mM)^{-1}','FontWeight','Bold','FontSize',12);
%set(gca,'YLim',[0 250])
return
