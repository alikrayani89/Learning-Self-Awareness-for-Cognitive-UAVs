function [outputArg1] = calculateLamdaS_basedOn_Maha(totNumOfSuperstates, measurement, ...
    meanOfSuperstates, R, covarianceOfSuperstates)
%% Input:
% totNumOfSuperstates: total number of Superstates
% measurement: the current observation
% meanOfSuperstates: matrix consisting of the mean value of each Superstate
% covarianceOfSuperstates: cell consisting of the covariance matrix of each Superstate

%% Calculate lambda in terms of battacharyya distance b/w observation & each superstate:
for index_s = 1:totNumOfSuperstates
    lamdaS(1, index_s) = MahalanobisDistance(measurement, ...
        meanOfSuperstates(index_s,:), topdm(R), topdm(covarianceOfSuperstates{1,index_s}));
end

%% Convert lamda to a discrete probability distribution:
n = 3; % using n can help to make the probability distribution more skewed (for example if n=1 give you [0.6 0.4], n=2 will give you [0.1 0.9])
probability_lamdaS = (1./((lamdaS).^n))/sum(1./((lamdaS).^n));

%% Output:
outputArg1 = probability_lamdaS;

end

