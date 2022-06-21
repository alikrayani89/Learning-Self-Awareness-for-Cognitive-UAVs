
% Input and Output: net structure obtained with GNG clustering
% This function calculates the overall transition matrix.
function [net] = CalculateTransitionMatrix (net)

%%  Transition matrix
transitionMat = zeros(net.N,net.N);
% Total length of training data
currLength = size(net.data, 1);
% Sequence of superstates of training data
nodesInTime = net.dataColorNode; 
trackDicrete = nodesInTime(1:currLength);
ind = find(diff(trackDicrete) ~= 0);

%%%%%
for k = 1:currLength-1
    transitionMat(trackDicrete(k,1),trackDicrete(k+1,1)) =...
        transitionMat(trackDicrete(k,1),trackDicrete(k+1,1)) + 1;
end
transitionMat = transitionMat./repmat(sum(transitionMat,2) + (sum(transitionMat,2)==0),1,net.N);

net.transitionMat = transitionMat;

end