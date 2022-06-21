function [net] = CalculateMaxClustersTime(net)

% Total length of training data
currLength = size(net.data, 1);
% Sequence of superstates of training data
zonesTrajectories = net.dataColorNode; 
% Number of clusters
N = net.N;

% Where to insert the max time values
maxClustersTime = zeros(1, N);

% Initialise the variable containing the previous zone of a comparison
% with the first zone of the signal
oldZone = zonesTrajectories(1,1);
finalZone = zonesTrajectories(currLength,1);
    
% Initialise the length of the run
currentMax = 0;
currentMax2 = 0;
for t = 1: currLength
    
    % Zone at current time
    newZone = zonesTrajectories(t);
    
    % If the zone has not changed
    if newZone ~= finalZone
        if newZone == oldZone
            % Increment the max of current run
            currentMax = currentMax + 1;
            % if the zone has changed
        else
            % %         oldZone = newZone;
            % Check if a longer run has been foung
            if currentMax > maxClustersTime(oldZone)
                maxClustersTime(oldZone) = currentMax;
            end
            oldZone = newZone;
            % Reinitialize current max
            currentMax = 1;
        end
    else
        if currentMax > maxClustersTime(oldZone)
                maxClustersTime(oldZone) = currentMax;
        end
        % Increment the max of current run
        currentMax2 = currentMax2 + 1;
        if currentMax > maxClustersTime(newZone)
            maxClustersTime(newZone) = currentMax2;
        end
    end
end

 

net.maxClustersTime = maxClustersTime;

end