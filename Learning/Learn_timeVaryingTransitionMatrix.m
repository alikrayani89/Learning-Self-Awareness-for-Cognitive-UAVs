function [outputArg1] = Learn_timeVaryingTransitionMatrix(netLTE, lengthInputData)
%LEARN_TIMEVARYINGTRANSITIONMATRIX Summary of this function goes here
%   Detailed explanation goes here

nodesInTime = netLTE.dataColorNode;
transitionMat = zeros(netLTE.N,netLTE.N);
timesSpent = [];
% % currLength = length(inputData);
currLength = lengthInputData;
trackDicrete = nodesInTime(1:currLength);
ind = find(diff(trackDicrete) ~= 0);

%% 3nd Method : Temporal Transition Matrix -> akr 27.03.2020 (FUNZIONA BENISSIMO)
codeInd = [0; ind];
diffe = size(netLTE.N,1)-codeInd(end,1);
if diffe>0
    codeInd = [codeInd; size(netLTE.N,1)];
end
tspentTran = diff(codeInd);
totCodeInd = size(codeInd,1);
kk = 2;
for k = 1:currLength
    timeMats{1, k} = zeros(netLTE.N, netLTE.N);
    timeMats{1,k}(trackDicrete(k,1), trackDicrete(k,1)) = 1;
    if kk <= totCodeInd
        if k == codeInd(kk,1)+1 % the next time instant the signal will transit to another super-state
            lastSuperState = trackDicrete(k-1,1);
            nextSuperState = trackDicrete(k,1);
            timeMats{1, k} = zeros(netLTE.N, netLTE.N);
            timeMats{1,k}(lastSuperState, lastSuperState) = tspentTran(kk-1,1);
            timeMats{1,k}(lastSuperState, nextSuperState) = 1;
            kk = kk + 1;
        end
    end
end

%% 4th Method : Temporal Transition Matrix (prob decrease/increase with time) -> akr 30.03.2020
% % codeInd = [0; ind];
% % diffe = size(netLTE.N,1)-codeInd(end,1);
% % if diffe>0
% %     codeInd = [codeInd; size(netLTE.N,1)];
% % end
% % tspentTran = diff(codeInd);
% % totCodeInd = size(codeInd,1);
% % kk = 1;
% % kkk = 0;
% % for k = 1:currLength
% %     timeMats{1, k} = zeros(netLTE.N, netLTE.N);
% %     timeMats{1,k}(trackDicrete(k,1), trackDicrete(k,1)) = tspentTran(kk, 1) - kkk;
% %     timeMats{1,k}(trackDicrete(k,1), trackDicrete(tspentTran(kk,1)+1, 1)) = kkk;
% %     kkk = kkk + 1;
% %         if kkk == tspentTran(kk, 1) % the next time instant the signal will transit to another super-state
% % %             lastSuperState = trackDicrete(k-1,1);
% % %             nextSuperState = trackDicrete(k,1);
% % %             timeMats{1, k} = zeros(netLTE.N, netLTE.N);
% % %             timeMats{1,k}(lastSuperState, lastSuperState) = tspentTran(kk-1,1);
% % %             timeMats{1,k}(lastSuperState, nextSuperState) = 1;
% %             kk = kk + 1;
% %             kkk = 0;
% %         end
% % end

netLTE.transMatsTimeAKR = timeMats;
outputArg1 = netLTE;

end

