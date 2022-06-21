function d=MahalanobisDistance(A, B, covA, covB)
% Return mahalanobis distance of two data matrices

% A is the mean of the first set (mean of x)
% B is the mean of the second set (mean of S)
% covA is the covariance of A (cov of x)
% covB is the covariance of B (cov of S)
n1 = 36;
n2 = 36;
n = n1 + n2;
xDiff=A-B;       % mean difference row vector

d=sqrt(xDiff*inv(covB)*xDiff'); % mahalanobis distance

end