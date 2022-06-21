function[tMax] = Find_Max_Time_Before_Transition(zonesSignal)

% Initialise the max to 1
tMax = 1;

% Initialise the variable containing the previous zone of a comparison
% with the first zone of the signal
oldZone = zonesSignal(1,1);
    
% Initialise the length of the run
currentMax = 1;
    
% Looping over the time instants of the trajectory
timeInstantsCurrentTrajectory = size(zonesSignal, 1);
for t = 2 : timeInstantsCurrentTrajectory

    % Zone at current time
    newZone = zonesSignal(t);

    % If the zone has not changed
    if newZone == oldZone
        % Increment the max of current run
        currentMax = currentMax + 1;
    % if the zone has changed
    else
        oldZone = newZone;
        % Check if a longer run has been foung
        if currentMax > tMax
            tMax = currentMax;
        end
        % Reinitialize current max
        currentMax = 1;
    end  
end

end