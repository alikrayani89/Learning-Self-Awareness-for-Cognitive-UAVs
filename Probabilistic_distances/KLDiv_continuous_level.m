function [outputArg1] = KLDiv_continuous_level(mu1, cov1, mu2, cov2)
%KLD_CONTINUOUS_LEVEL Summary of this function goes here
%   Detailed explanation goes here

cov1 = topdm(cov1); % topdm check and then convert
cov2 = topdm(cov2); % topdm check and then convert

%% Dkl = 1/2(tr(cov2^(-1)cov1)+(mu2-mu1)^T)
k = size(mu1,1); % data dimensionality

%% KLD(p||q)
% % Dkl = (1)*(trace(inv(cov2).*cov1) + transpose(mu2-mu1)*inv(cov2)*(mu2-mu1) - k + log(det(cov2)./det(cov1)));

%% KLD(q||p)
% % Dkl = (1)*(trace(inv(cov1).*cov2) + transpose(mu1-mu2)*inv(cov1)*(mu1-mu2) - k + log(det(cov1)./det(cov2)));

%% Symmetric KLD
%% (1)*KLD(p||q) + (1)*KLD(q||p) or (1/2)*KLD(p||q) + (1/2)*KLD(q||p)
Dkl = (1/2)*(trace(inv(cov2).*cov1) + transpose(mu2-mu1)*inv(cov2)*(mu2-mu1) - k + log(det(cov2)./det(cov1)))+...
    (1/2)*(trace(inv(cov1).*cov2) + transpose(mu1-mu2)*inv(cov1)*(mu1-mu2) - k + log(det(cov1)./det(cov2)));

outputArg1 = Dkl;
end

