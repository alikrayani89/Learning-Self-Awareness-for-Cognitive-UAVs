function d=bhattacharyyadistance(mu1, mu2, C1, C2)
% BHATTACHARYYA  Bhattacharyya distance between two Gaussian classes
% Inputs: mu1 mu2 mean c1 c2 covariance
% Reference:
% Kailath, T., The Divergence and Bhattacharyya Distance Measures in Signal
% Selection, IEEE Trasnactions on Communication Technology, Vol. 15, No. 1,
% pp. 52-60, 1967
%
% By Yi Cao at Cranfield University on 8th Feb 2008.
%
C=(C1+C2)/2;
dmu=(mu1-mu2)/chol(C);
% % dmu=(mu1-mu2)/lu(C); % akr: or qr ma chol molto meglio
try
    d=0.125*dmu*dmu'+0.5*log(det(C/chol(C1*C2)));
catch
    d=0.125*dmu*dmu'+0.5*log(abs(det(C/sqrtm(C1*C2))));
    %warning('MATLAB:divideByZero','Data are almost linear dependent. The results may not be accurate.');
end
% d=0.125*dmu*dmu'+0.25*log(det((C1+C2)/2)^2/(det(C1)*det(C2)));