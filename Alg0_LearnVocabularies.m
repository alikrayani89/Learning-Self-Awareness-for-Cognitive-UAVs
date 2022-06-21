% Function to perform clustering

clc
clear
close all

curDir = pwd;
cd(curDir)

addpath('./Learning');
addpath('./Plot_scripts');

%% Select a certain example from DataSet
example = 1;
snr = 20;

dirtemp = strcat(curDir, '\DataSet\Ex', num2str(example), '\SNR_', num2str(snr), 'dB');
cd(dirtemp)

%% Load signal
load('DataX_LTE')
training_states = DataX_LTE;

output_GNG = 'Vocabulary_.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading of training data
training_GS = [training_states];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perfoming clustering with the GNG algorithm
%% GNG parameters 1
params.N = 10;                                                               %    Number of nodes
params.MaxIt = 50;                                                           %    Iteration (repetition of input data)
params.L_growing = 1000;                                                     %    Growing rate
params.epsilon_b = 0.05;                                                     %    Movement of winner node
params.epsilon_n = 0.0006;                                                   %    Movement of all other nodes except winner
params.alpha = 0.5;                                                          %    Decaying global error and utility
params.delta = 0.9995;                                                       %    Decaying local error and utility
params.T = 100;                                                              %    It could be a function of params.L_growing, e.g., params.LDecay = 2*params.L_growing
params.L_decay = 1000;                                                       %    Decay rate sould be faster than the growing then it will remove extra nodes
params.alpha_utility = 0.0005;                                               %    It could be a function of params.delta, e.g., params.alpha_utility = 0.9*params.delta
params.k = 10;
params.seedvector = 1;

% Input to GNG == GS states
inputData = training_GS;

%% Calling GNG
% % net = GrowingNeuralGasNetwork(InputData, params, false);
plotflagGNG = true;
GSVdimensions = size(training_GS, 2);
net = GrowingNeuralGasNetwork(inputData, params, plotflagGNG, GSVdimensions);

%% Plot >:
plotClusters(inputData, net, '')
% % plotMap = true;
%% Plot <:

%% Calculate transition matrix
net = CalculateTransitionMatrix(net);

%% Plot transition matrix
figure
imagesc(net.transitionMat);
colorbar;
colormap jet
title (strcat('Transition Matrix ', ''), 'interpreter','latex')

%% Plot transition Matrix as a Graph
mc = dtmc(net.transitionMat);
figure;
graphplot(mc,'ColorEdges',true);
title (strcat('Graph Transition Matrix ', ''), 'interpreter','latex')


%% Calculate temporal transition matrix
net = CalculateTemporalTransitionMatrices (net);

net = Learn_timeVaryingTransitionMatrix(net, length(training_GS));

%% Find the max time spent in each cluster
net = CalculateMaxClustersTime (net);

%% Learn Conditional (transitional) Generalized Mean Values (GMC)
net = DiscoverGeneralizedMeanValues (net);

%% Plot GMCs to evaluate the output:
plotGeneralizedMeanValues(net);

% Saving the GNG output
save(output_GNG, 'net');

