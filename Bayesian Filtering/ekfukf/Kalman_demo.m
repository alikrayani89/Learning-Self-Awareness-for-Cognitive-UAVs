clear
% Stepsize
dt = 0.5;

% Transition matrix for the continous-time system.
A = [1  0  0  0;
     0  1  0  0;
     0  0  1  0;
     0  0  0  1];

% Process noise variance
sig_a = 0.2;
Q = eye(4) * sig_a^2;

% Measurement model.
H = [1 0 0 0;
     0 1 0 0];

% Variance in the measurements.
r1 = 1;
R = diag([r1 r1]);



M = 100; % number of Monte-Carlo runs

for mnt=1:M
    
data_generator % generate a trajectory
load('data.mat');

clear MM PP m P

% Initial guesses for the state mean and covariance.
m = [0 0 0 0]';
P = diag([0.1 0.1 0.1 0.1]);

% Space for the estimates.
MM = zeros(size(m,1), size(Y,2));
PP = zeros(size(m,1), size(m,1), size(Y,2));


% Filtering steps.
for i = 1:size(Y,2)
   [m,P] = kf_predict(m,P,A,Q);
   [m,P] = kf_update(m,P,Y(:,i),H,R);
   MM(:,i) = m;
   PP(:,:,i) = P;
end

[mse, X_r] = compare(MM); % compare with ground-truth
all_mse(:,mnt) = mse; % save MSE values
% mnt


end

mean(all_mse,2) % average MSE values

% plot last run
subplot(1,2,1);
plot(X_r(1,:), X_r(2,:),'--', MM(1,:), MM(2,:),X_r(1,1),X_r(2,1),...
     'o','MarkerSize',12)
legend('Real trajectory', 'Filtered');
title('Position estimation with Kalman filter.');
xlabel('x');
ylabel('y');

subplot(1,2,2);
plot(X_r(3,:), X_r(4,:),'--', MM(3,:), MM(4,:),X_r(3,1),...
     X_r(4,1),'ro','MarkerSize',12);
legend('Real velocity', 'Filtered');
title('Velocity estimation with Kalman filter.');
xlabel('x^.');
ylabel('y^.');