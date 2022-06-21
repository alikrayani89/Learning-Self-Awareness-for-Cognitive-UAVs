
% Input and Output: net structure obtained with GNG clustering
% The function creates a set of transition matrices, one for each time value
% from t to tMax, being tMax the time we have already spent in a node. 
function [net] = CalculateTemporalTransitionMatrices (net)

%% Temporal Transition Matrix

% Total length of training data
currLength = size(net.data, 1);
% Sequence of superstates of training data
nodesInTime = net.dataColorNode; 
% Number of clusters
N = net.N;

% Find max number of time instants before a zone change 
tMax = Find_Max_Time_Before_Transition(nodesInTime);

% Find TIME TRANSITION MATRICES for follower
transMatsTime = Find_Time_Transition_Matrices(tMax, ...
    N, nodesInTime, currLength);

net.transMatsTime = transMatsTime;

end