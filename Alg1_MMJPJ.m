
clc
clear
close all

curDir = pwd;
cd(curDir)

addpath('./Testing_codes');
addpath('./Probabilistic_distances')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters

%% Select a certain example from DataSet
example = 1;
snr = 18;

dirtemp = strcat(curDir, '\DataSet\Ex', num2str(example), '\SNR_', num2str(snr), 'dB');
cd(dirtemp)

% % load('DataX_LTE_Jammed_BPSK')
% % inputData = DataX_LTE_Jammed_BPSK;

% % load('DataX_LTE_Jammed_8PSK')
% % inputData = DataX_LTE_Jammed_8PSK;

load('DataX_LTE_Jammed_256QAM')
inputData = DataX_LTE_Jammed_256QAM;

vocabularyName = 'Vocabulary_.mat';

% Loading the vocabulary
load(vocabularyName);
Vocabulary = net;

dataTest = inputData';

N = 10; % total number of particles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALLING THE MJPF
[estimationAbn1] = M_MJPF(dataTest, Vocabulary, N, curDir);

