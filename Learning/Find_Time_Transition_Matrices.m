function[transMatsTime] = Find_Time_Transition_Matrices(tMax, somSize, zonesSignal, ...
    numberTrainingData)


transMatsTime = cell(1, tMax);
for i = 1: tMax
    transMatsTime{1, i} = zeros(somSize);
end


% I take the first zone
prevZone = zonesSignal(1, 1);
time = 0;

% loop over the number of time instants of the data signal
for t=2:numberTrainingData
    % I add a time instant
    time = time + 1;
    if (time > tMax)
        break
    end
    % New zone
    newZone = zonesSignal(t);
    % If I change the zone with respect to the previous time
    % instant(/s)
    if (prevZone ~= newZone)
        % I increment the corresponding value in the correct transition
        % matrix
        transMatsTime{1,time}(prevZone, newZone) = ...
            transMatsTime{1,time}(prevZone, newZone) + 1;
        % And I update the zone value
        prevZone = newZone;

        % I reinitialize the time
        time = 0;
    % Otherwise, if I remain in the same zone
    else
        transMatsTime{1,time}(prevZone, prevZone) = ...
            transMatsTime{1,time}(prevZone, prevZone) + 1;
    end

end

% For each transition matrix
for t = 1:tMax
    % looping over rows of current matrix
    for row= 1: somSize
        % sum of the elements of the row
        sumElementsRow = sum(transMatsTime{1, t}(row, :));
        
        % to prevent division by 0
        if sumElementsRow ~=0
            % looping over columns of current matrix
            for column = 1: somSize
                % normalise matrix element
                transMatsTime{1,t}(row, column) = ...
                    transMatsTime{1,t}(row, column)/sumElementsRow;
            end
        end
    end
end

end