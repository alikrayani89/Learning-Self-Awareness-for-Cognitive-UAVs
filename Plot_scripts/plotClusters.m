function [] = plotClusters(inputData, net, stringModulationFormat)
%PLOTCLUSTERS Summary of this function goes here
%   Detailed explanation goes here

mycolors = colorcube;
mycolors = [mycolors;mycolors];

h = figure;
hold on
scatter(inputData(:,1),inputData(:,2),60,mycolors(net.dataColorNode,:),'.')
scatter(net.w(:,1),net.w(:,2),250,'+','k','linewidth',2)
grid on
box on
xlim([-1.5 1.5])
ylim([-1.5 1.5])
title (stringModulationFormat, 'interpreter','latex')
% % xlabel('In-Phase (I)', 'FontSize',14, 'interpreter','latex')
% % ylabel('Quadrature (Q)', 'FontSize',14, 'interpreter','latex')
xlabel('$\bf{In-Phase \ (I)}$', 'FontSize',12, 'interpreter','latex');
ylabel('\bf{Quadrature (Q)}', 'FontSize',12, 'interpreter','latex');
%     saveas(gcf, 'f132QAM', 'bmp')

h1 = figure;
hold on
scatter(inputData(:,1),inputData(:,2),60,mycolors(net.dataColorNode,:),'.')
scatter(net.w(:,1),net.w(:,2),200'.','k','linewidth',2)
quiver(net.w(:,1), net.w(:,2), net.w(:,3), net.w(:,4),'k','Autoscale','off','LineWidth',1)
grid on
box on
xlim([-1.5 1.5])
ylim([-1.5 1.5])
title (stringModulationFormat, 'interpreter','latex')
xlabel('$\bf{In-Phase \ (I)}$', 'FontSize',12, 'interpreter','latex');
ylabel('\bf{Quadrature (Q)}', 'FontSize',12, 'interpreter','latex');
%     saveas(gcf, 'f232QAM', 'bmp')

end
