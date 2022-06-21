function [outputArg1] = KLD_Abnormality(totNumOfSuperstates, N, histogram,...
    transitionMat, probability_lamdaS)

%% Input:
% totNumOfSuperstates: total number of Superstates
% N: total number of Particles
% histogram: histogram at time t-1 (after PF resampling) 
% transitionMat: the transition matrix learned from previous experience
% probability_lamdaS: probability vector representing a discrete probability disctribution

%% Procedure:
sommaKLD_simmetrica = 0;
for indKLD = 1:totNumOfSuperstates
    particella = histogram(1,indKLD);
    if particella>0
        PP = transitionMat(indKLD,:)+1e-100; % add 1e-100 since KLD doesnt allow zero values
        QQ = probability_lamdaS;
        KLD_simmetrica = (particella/N)*KLDiv(PP,QQ) + (particella/N)*KLDiv(QQ,PP); %to achieve symmerty
        sommaKLD_simmetrica = sommaKLD_simmetrica + KLD_simmetrica;
    end
end

%% Output:
outputArg1 = sommaKLD_simmetrica;

end

