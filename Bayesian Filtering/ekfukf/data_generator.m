function data_generator

% Stepsize
dt = 0.5;

% $x_k = Ax_{k-1}+Ga_k$
% Transition matrix for the continous-time system.
A = [1  0 dt  0;
     0  1  0 dt;
     0  0  1  0;
     0  0  0  1];
% acceleration term
G = [dt^2/2 0 dt 0;
     0 dt^2/2 0 dt]';

% Process noise variance
sig_a = 0.2;
Qa = eye(2)*sig_a^2;
Q = G*Qa*G';

% Measurement model.
H = [1 0 0 0;
     0 1 0 0];

% Variance in the measurements.
r1 = 10;
R = diag([r1 r1]);

% Generate the data.
n = 500;
Y = zeros(size(H,1),n);
X_r = zeros(size(A,1),n);
X_r(:,1) = [0 0 1 1]';
for i = 2:n
   X_r(:,i) = A*X_r(:,i-1) + G*gauss_rnd(zeros(size(Qa,1),1), Qa);
end

% Generate the measurements.
for i = 1:n
   Y(:,i) = H*X_r(:,i) + gauss_rnd(zeros(size(Y,1),1), R);
end

plot(X_r(1,:),X_r(2,:),Y(1,:),Y(2,:),'.',X_r(1,1),...
     X_r(2,1),'ro','MarkerSize',12);

save('data.mat', 'Y')
save('gt.mat', 'X_r')