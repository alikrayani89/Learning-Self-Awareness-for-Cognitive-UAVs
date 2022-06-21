function [mse, X_r] = compare(Xf)
load('gt.mat')

mse = mean((Xf-X_r).^2,2);