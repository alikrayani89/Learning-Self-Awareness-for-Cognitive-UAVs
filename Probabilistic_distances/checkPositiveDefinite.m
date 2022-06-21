function [outputArg1] = checkPositiveDefinite(covariance)
%CHECKPOSITIVEDEFINITE Summary of this function goes here
%   Detailed explanation goes here

%% A matrix is positive definite if all it's associated eigenvalues are positive.
%% A way to check if a certain matrix is positive definite is as follows:

eig_covariance = eig(covariance);
flag = 0;

for i = 1:rank(covariance)
    if eig_covariance(i) <= 0
        flag = 1;
    end
end

% % if flag == 1
% %   disp('the matrix is not positive definite')
% %   else
% %   disp('the matrix is positive definite')
% % end

outputArg1 = flag;

end

