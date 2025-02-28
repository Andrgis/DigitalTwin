%% Unified Hybrid Fault Diagnosis Simulation for an ASV
% This script simulates the dynamics of a small ASV (the Otter) under sensor and actuator faults.
% It implements a simplified adaptive extended Kalman filter (AEKF) for joint state and fault estimation.
% The simulation parameters, state equations, sensor filter equations, and fault injection times are chosen to
% mimic the example presented in the paper.
%
% References:
% A. Hasan and P.S. Rossi, "A unified sensor and actuator fault diagnosis in digital twins for remote operations,"
% Mechanical Systems and Signal Processing, vol. 222, 2025. :contentReference[oaicite:0]{index=0}

clear; clc; close all;

%% Simulation Parameters
T = 45;              % total simulation time in seconds
dt = 0.01;           % sampling time (s)
N = round(T/dt);     % number of simulation steps
time = (dt:dt:T)';    % time vector

% Filter tuning parameters (example values)
mu = 0.5;            % filtering matrix parameter for sensor filter
nu = 0.5;            % scaling parameter in sensor fault injection
lambda = 0.95;       % forgetting factor for fault adaptation
gamma = 0.01;        % learning rate for fault estimation update

%% True System and Fault Initialization
% True state: x = [p; q; U; chi]
x_true = [0; 0; 1; 0];  % initial state (p=0, q=0, U=1, chi=0)

% Preallocate arrays for storing histories
x_true_hist    = zeros(4, N);
x_est_hist     = zeros(4, N);
fault_est_hist = zeros(6, N);  % fault estimates: [theta_a1; theta_a2; theta_s1; theta_s2; theta_s3; theta_s4]
z_meas_hist    = zeros(4, N);
u_hist         = zeros(2, N);  % control inputs: [a; r]

% Initialize true fault signals (piecewise constant)
% Actuator faults: theta_a1 and theta_a2 (faults reduce effectiveness)
theta_a1_true = zeros(N,1);
theta_a2_true = zeros(N,1);
% Sensor faults: theta_s1, theta_s2, theta_s3, theta_s4 (bias faults)
theta_s1_true = zeros(N,1);
theta_s2_true = zeros(N,1);
theta_s3_true = zeros(N,1);
theta_s4_true = zeros(N,1);

for k = 1:N
    t = (k-1)*dt;
    % Actuator faults: at t >= 10 s, 10% reduction; at t >= 20 s, 20% reduction.
    if t >= 10 && t < 20
        theta_a1_true(k) = 0.1;
        theta_a2_true(k) = 0.1;
    elseif t >= 20
        theta_a1_true(k) = 0.2;
        theta_a2_true(k) = 0.2;
    end
    % Sensor faults: at t >= 20 s, introduce a bias of 0.05; at t >= 30 s, increase bias to 0.1.
    if t >= 20 && t < 30
        theta_s1_true(k) = 0.05;
        theta_s2_true(k) = 0.05;
        theta_s3_true(k) = 0.05;
        theta_s4_true(k) = 0.05;
    elseif t >= 30
        theta_s1_true(k) = 0.1;
        theta_s2_true(k) = 0.1;
        theta_s3_true(k) = 0.1;
        theta_s4_true(k) = 0.1;
    end
end

% Noise covariances
Q_proc = diag([0.01, 0.01, 0.01, 0.01]);  % process noise covariance (state)
R_meas = diag([0.05, 0.05, 0.05, 0.05]);  % measurement noise covariance (sensor)

%% Filter Initialization
% Initial state estimate (assume perfect knowledge initially)
x_est = x_true;
% Initial fault estimates (assume no fault initially)
theta_est = zeros(6,1);   % [theta_a1; theta_a2; theta_s1; theta_s2; theta_s3; theta_s4]
% Initial state covariance
P = eye(4)*0.1;
% Initialize sensor filter output (used to mimic sensor signal processing)
z_est = x_true;  % start with the true state as the initial sensor filter output

% Save initial histories
x_true_hist(:,1)    = x_true;
x_est_hist(:,1)     = x_est;
fault_est_hist(:,1) = theta_est;
z_meas_hist(:,1)    = z_est;

%% Main Simulation Loop
for k = 2:N
    t = (k-1)*dt;
    
    %% Control Inputs
    % For illustration, use constant control inputs:
    a = 1;      % acceleration (m/s^2)
    r = 0.1;    % turning rate (rad/s)
    u = [a; r];
    u_hist(:,k) = u;
    
    %% True System Dynamics
    % Retrieve current actuator fault values
    theta_a1 = theta_a1_true(k-1);
    theta_a2 = theta_a2_true(k-1);
    
    % Generate process noise (sampled from a Gaussian distribution)
    w = sqrt(diag(Q_proc)).*randn(4,1);
    
    % Update state using the nonlinear model:
    % p(k) = p(k-1) + dt * U*cos(chi)
    % q(k) = q(k-1) + dt * U*sin(chi)
    % U(k) = U(k-1) + dt*(1 - theta_a1)*a
    % chi(k) = chi(k-1) + dt*(1 - theta_a2)*r
    p   = x_true(1) + dt*x_true(3)*cos(x_true(4)) + dt*w(1);
    q   = x_true(2) + dt*x_true(3)*sin(x_true(4)) + dt*w(2);
    U_  = x_true(3) + dt*(1 - theta_a1)*a + dt*w(3);
    chi = x_true(4) + dt*(1 - theta_a2)*r + dt*w(4);
    x_true = [p; q; U_; chi];
    x_true_hist(:,k) = x_true;
    
    %% Sensor Measurements via a Filtering Process
    % Retrieve current sensor fault values
    theta_s1 = theta_s1_true(k-1);
    theta_s2 = theta_s2_true(k-1);
    theta_s3 = theta_s3_true(k-1);
    theta_s4 = theta_s4_true(k-1);
    
    % Generate measurement noise
    v = sqrt(diag(R_meas)).*randn(4,1);
    
    % The sensor filter (mimicking the paper’s approach):
    % z(k) = mu*dt*x_true(k-1) + (1 - mu*dt)*z(k-1) + mu*nu*dt*theta_s + mu*dt*v
    z = mu*dt*x_true_hist(1:4,k-1) + (1 - mu*dt)*z_est + mu*nu*dt*[theta_s1; theta_s2; theta_s3; theta_s4] + mu*dt*v;
    z_est = z;  % update sensor filter output
    z_meas_hist(:,k) = z;
    
    %% Adaptive Extended Kalman Filter (AEKF) for State Estimation
    % The state prediction uses the current fault estimates (for actuator faults)
    theta_a_est = theta_est(1:2);
    x_pred = zeros(4,1);
    x_pred(1) = x_est(1) + dt*x_est(3)*cos(x_est(4));
    x_pred(2) = x_est(2) + dt*x_est(3)*sin(x_est(4));
    x_pred(3) = x_est(3) + dt*(1 - theta_a_est(1))*a;
    x_pred(4) = x_est(4) + dt*(1 - theta_a_est(2))*r;
    
    % Compute the Jacobian F = df/dx evaluated at the current estimate
    F = eye(4);
    F(1,3) = dt*cos(x_est(4));
    F(1,4) = -dt*x_est(3)*sin(x_est(4));
    F(2,3) = dt*sin(x_est(4));
    F(2,4) = dt*x_est(3)*cos(x_est(4));
    % F(3,:) and F(4,:) remain as in the linearized update.
    
    % State covariance prediction
    P = F * P * F' + Q_proc;
    
    % Predicted measurement using the sensor filter model
    % h(x, theta_s) = mu*dt*x_pred + (1 - mu*dt)*z_est + mu*nu*dt*theta_s_est
    theta_s_est = theta_est(3:6);
    z_pred = mu*dt*x_pred + (1 - mu*dt)*z_est + mu*nu*dt*theta_s_est;
    
    % Measurement Jacobian (approximate as H = mu*dt*I)
    H = mu*dt*eye(4);
    
    % Innovation (measurement residual)
    y_tilde = z - z_pred;
    
    % Innovation covariance
    S_meas = H * P * H' + R_meas;
    
    % Kalman gain for state update
    K = P * H' / S_meas;
    
    % Update state estimate and covariance
    x_est = x_pred + K * y_tilde;
    P = (eye(4) - K * H) * P;
    
    %% Fault Estimation Update (Dual Estimation)
    % Here we use a simple gradient update for the fault estimates based on the observed residuals.
    % For actuator faults: use the error in the state prediction of U and chi.
    fault_error_a = [ (x_est(3) - x_pred(3)) / (-dt*a);
                      (x_est(4) - x_pred(4)) / (-dt*r) ];
    theta_a_est = theta_a_est + gamma * fault_error_a;
    
    % For sensor faults: update based on the measurement residual
    fault_error_s = y_tilde / (mu*nu*dt);
    theta_s_est = theta_s_est + gamma * fault_error_s;
    
    % Combine fault estimates
    theta_est = [theta_a_est; theta_s_est];
    fault_est_hist(:,k) = theta_est;
    
    % Save updated state estimate
    x_est_hist(:,k) = x_est;
end

%% Plot Results

% Plot true vs. estimated state variables
figure;
subplot(2,2,1);
plot(time, x_true_hist(1,:), 'b', time, x_est_hist(1,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('p');
legend('True','Estimated'); title('State: p');

subplot(2,2,2);
plot(time, x_true_hist(2,:), 'b', time, x_est_hist(2,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('q');
legend('True','Estimated'); title('State: q');

subplot(2,2,3);
plot(time, x_true_hist(3,:), 'b', time, x_est_hist(3,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('U');
legend('True','Estimated'); title('State: U');

subplot(2,2,4);
plot(time, x_true_hist(4,:), 'b', time, x_est_hist(4,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('\chi');
legend('True','Estimated'); title('State: \chi');

% Plot actuator fault estimates
figure;
subplot(2,1,1);
plot(time, theta_a1_true, 'b', time, fault_est_hist(1,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('\theta_{a1}');
legend('True','Estimated'); title('Actuator Fault: \theta_{a1}');

subplot(2,1,2);
plot(time, theta_a2_true, 'b', time, fault_est_hist(2,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('\theta_{a2}');
legend('True','Estimated'); title('Actuator Fault: \theta_{a2}');

% Plot sensor fault estimates
figure;
subplot(2,2,1);
plot(time, theta_s1_true, 'b', time, fault_est_hist(3,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('\theta_{s1}');
legend('True','Estimated'); title('Sensor Fault: \theta_{s1}');

subplot(2,2,2);
plot(time, theta_s2_true, 'b', time, fault_est_hist(4,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('\theta_{s2}');
legend('True','Estimated'); title('Sensor Fault: \theta_{s2}');

subplot(2,2,3);
plot(time, theta_s3_true, 'b', time, fault_est_hist(5,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('\theta_{s3}');
legend('True','Estimated'); title('Sensor Fault: \theta_{s3}');

subplot(2,2,4);
plot(time, theta_s4_true, 'b', time, fault_est_hist(6,:), 'r--', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('\theta_{s4}');
legend('True','Estimated'); title('Sensor Fault: \theta_{s4}');
