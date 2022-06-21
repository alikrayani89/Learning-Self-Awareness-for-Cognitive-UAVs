function [netLTE] = DiscoverGeneralizedMeanValues(netLTE)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Mean of transitions between node: START >
nodesMean_Transitions = cell(1, netLTE.N);
nodesCov_Transitions = cell(netLTE.N, netLTE.N);

num_notDerivativeElements = size(netLTE.data,2)/2;
num_DerivativeElements = num_notDerivativeElements;
num_totalElements = num_notDerivativeElements + num_DerivativeElements;
for j = 1:netLTE.N
    for k = 1:netLTE.N
        if j==k
            nodesMean_Transitions{1,j}(k,:) =  netLTE.nodesMean(k,:);
            nodesCov_Transitions{j,k} = netLTE.nodesCov{1, k};
        else
            for z=1:1:num_notDerivativeElements
                mu(1,z) = netLTE.nodesMean(k,z) - netLTE.nodesMean(j,z);
                nodesMean_Transitions{1,j}(k,z) = 0;
            end
            indexZ = 1;
            for L = num_DerivativeElements+1:1:num_totalElements
                nodesMean_Transitions{1,j}(k,L) = mu(1,indexZ);
                indexZ = indexZ + 1;
            end
            clear indexZ;
            nodesCov_Transitions{j,k} = netLTE.nodesCov{1, j} + netLTE.nodesCov{1, k};
        end
    end
end

netLTE.nodesMean_Transitions = nodesMean_Transitions;
netLTE.nodesCov_Transitions = nodesCov_Transitions;
%% Mean of transitions between node: END >

%% NewMean of each cluster >
for i = 1:netLTE.N
    Iclean_signal = netLTE.datanodes{1,i}(:,1:num_notDerivativeElements/2);
    Qclean_signal = netLTE.datanodes{1,i}(:,(num_notDerivativeElements/2)+1:(num_totalElements/2));
    
    total_ofdm_symbols = size(netLTE.datanodes{1,i}, 1);
    total_frequencies = num_DerivativeElements/2;
    % I components
    IdotClean = zeros(total_ofdm_symbols, total_frequencies);
    IdotClean_diff = diff(Iclean_signal);
    IdotClean(2:end, :) = IdotClean_diff;
    % Q components
    QdotClean = zeros(total_ofdm_symbols, total_frequencies);
    QdotClean_diff = diff(Qclean_signal);
    QdotClean(2:end, :) = QdotClean_diff;
    DataX_GS = [Iclean_signal, Qclean_signal, IdotClean, QdotClean];
    netLTE.datanodesNew{1,i} = DataX_GS;
    netLTE.nodesMeanNew(i,:) = mean(netLTE.datanodesNew{1,i});
    netLTE.nodesCovNew{1,i} = cov(netLTE.datanodesNew{1,i});
end
%% NewMean of each cluster >

end

