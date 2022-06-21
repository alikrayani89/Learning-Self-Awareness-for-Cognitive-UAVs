function [] = plotGeneralizedMeanValues(net)
%PLOTGENERALIZEDMEANVALUES Summary of this function goes here
%   Detailed explanation goes here

% 1) first plot the mean value of each cluster
%
for i = 1:1:net.N
    % % scatterplot(complex(net.nodesMean(i,1), net.nodesMean(i,2)), marker, 'x', color, 'b','linewidths', 1);
    scatter(net.nodesMean(i,1), net.nodesMean(i,2),'+','k','LineWidth',1.5);
    % % str = '$S_{ii}^{g1}$';
    str = strcat('$\tilde{S}_{', num2str(i), '}^{g_{1}}$');
    text(net.nodesMean(i,1)+0.04, net.nodesMean(i,2)+0.04, str, 'FontSize',14, 'interpreter','latex')
    % % text(net.nodesMean(i,1), net.nodesMean(i,2), str)
    box on
    xlabel('In-Phase (I)')
    ylabel('Quadrature (Q)')
    hold on
    % % i_prev = i;
    % 2) then plot the force learned if we stay in the same cluster mu_{ii} or-
    % if we pass to a new cluster mu_{ij}.
    % User quiver(x, y, \dot{x}, \dot{y}), where dotx and doty are the force or velocity
    for j = 1:1:net.N
        strMu = strcat('$\tilde{\mu}_{', num2str(i), num2str(j), '}^{g_{1}}$');
        if i == j
            quiver(net.nodesMean(i,1), net.nodesMean(i,2), ...
                net.nodesMean_Transitions{1, i}(j,3), net.nodesMean_Transitions{1, i}(j,4), 0.1);
        else
            quiver(net.nodesMean(i,1), net.nodesMean(i,2), ...
                net.nodesMean_Transitions{1, i}(j,3), net.nodesMean_Transitions{1, i}(j,4));
        end
        text(net.nodesMean_Transitions{1, i}(j,3), net.nodesMean_Transitions{1, i}(j,4), strMu, 'FontSize',14, 'interpreter','latex')
    end
end

% % % 2) then plot the force learned if we stay in the same cluster mu_{ii} or-
% % % if we pass to a new cluster mu_{ij}.
% % % User quiver(x, y, \dot{x}, \dot{y}), where dotx and doty are the force or velocity
% % quiver(net.nodesMean(i,1), net.nodesMean(i,2), mu12(i,1), mu12(i,2))

end

