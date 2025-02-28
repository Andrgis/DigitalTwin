%% Otter Actuator Fault Diagnosis using an Adaptive Extended Kalman Filter (AEKF)
clc; clear; clf;

%% Time horizon and simulation parameters
dt  = 0.001;
T   = 4;              % total simulation time (s)
t   = dt:dt:T;
N   = length(t);

%% Parameters from Table 1
Xu = -0.7225; Xuu = -1.3274; m = 23.8; X_udot = -2.0;
Yv = -0.8612; Yvv = -36.2823; Yv_dot = -10.0; xg = 0.046;
Yr = 0.1079; Nrr = -1; Yr_dot = 0.0; Nv_dot = 0.0;
Nv = 0.1052; Nr = -0.5; Izz = 1.76; Nr_dot = -1.0;

%% Initial state and fault parameters
% State: x = [eta; nu] = [x; y; yaw; u; v; r]
x     = [0; 0; 0; 100; 0.01; 0];  % true initial state
xhat  = x;                        % initial state estimate

% Actuator (fault) parameters: true fault is zero initially
theta     = zeros(3,1);           % true fault vector (unknown)
thetahat  = zeros(3,1);           % fault estimate

%% Noise covariances (tuning parameters)
n_x = length(x);
n_y = n_x;   % assume full-state measurements
QF = 0.01*eye(n_x);   % process noise covariance
RF = 0.04*eye(n_y);   % measurement noise covariance

%% AEKF initialization
Pplus       = eye(n_x);    % state error covariance
S           = 0.1*eye(3);  % fault-related covariance matrix
UpsilonPlus = zeros(n_x,3);
lambda      = 0.995;
a           = 0.999;

%% Matrices

A = eye(6);
% Mass matrix, Coriolis, and damping (as in your simulation code)
M = [ m - X_udot,           0,            0;
          0,     m - Yv_dot,  m*xg - Yr_dot;
          0,     m*xg - Nv_dot, Izz - Nr_dot ];
      
m11 = M(1,1);
m22 = M(2,2);
m23 = M(2,3);
m32 = M(3,2);
m33 = M(3,3);

%% Control input
% Nominal control input (for example, constant control input)
B_c = [zeros(3); 
       eye(3)/M];
tau_c = [-300; 100; 0 ];

% We choose Ψ = -B_tau; here B_tau = eye(3) in the dynamics affecting ν.
Psi = -eye(3)*diag(tau_c);

%% Storage for plotting
xArray      = zeros(n_x, N);
xhatArray   = zeros(n_x, N);
thetaArray  = zeros(3, N);
thetahatArray = zeros(3, N);
uArray      = zeros(3, N);

%% Simulation loop for fault diagnosis
for k = 1:N
    
    % ----- True System Dynamics -----
    % Compute the rotation matrix based on current yaw
    yaw = x(3);
    R = [cos(yaw), -sin(yaw), 0;
         sin(yaw),  cos(yaw), 0;
              0,         0, 1];
    
    nu = x(4:6);
    C_M = [ 0, 0, -M(2,2)*nu(2) - 0.5*(M(2,3)+M(3,2))*nu(3);
            0, 0, M(1,1)*nu(1);
            M(2,2)*nu(2) + 0.5*(M(2,3)+M(3,2))*nu(3), -M(1,1)*nu(1), 0 ];
    D_M = -[ Xu+Xuu*abs(nu(1)),       0,            0;
               0,       Yv+Yvv*abs(nu(2)),       Yr;
               0,               Nv,    Nr+Nrr*abs(nu(3)) ];
    
    % Here, we include the fault term in the control input:
    % tau_actual = tau_c + Ψ*theta.
    % (For simulation purposes, we can introduce a fault at a given time.)
    if k*dt > 2  % introduce fault after 2 seconds
        theta = [0.05; -0.03; 0.02];  % true fault (for example)
    end
    u_actual = tau_c + Psi*theta;
    
    % Update state using Euler integration for the Otter dynamics:
    % η̇ = R*ν,  and  ν̇ = inv(M)* ( - (C_M+D_M)*ν + B_tau*u_actual )
    eta_dot = R * nu;
    nu_dot  = inv(M)* ( - (C_M + D_M)*nu + u_actual );
    x(1:3) = x(1:3) + dt * eta_dot;
    x(4:6) = x(4:6) + dt * nu_dot;
    
    % Add process noise (optional)
    x = x + sqrt(dt)*mvnrnd(zeros(n_x,1), QF)';
    
    % ----- Measurement -----
    y = x + sqrt(dt)*mvnrnd(zeros(n_x,1), RF)';  % full state measurement
    
    % ----- Adaptive Extended Kalman Filter (AEKF) Update -----
    % Compute Jacobian of the system dynamics with respect to state (FX)
    % Here, we linearize the dynamics around the current estimate xhat.
    % For brevity, we approximate FX by perturbing the state or using an analytical Jacobian.
    % In this example, we assume a simplified FX (identity plus a dt term).
    FX = A+dt*[0 0 -sin(xhat(3))*xhat(4)-cos(xhat(3))*xhat(5) cos(xhat(3)) -sin(xhat(3)) 0;
        0 0 cos(xhat(3))*xhat(4)-sin(xhat(3))*xhat(5) sin(xhat(3)) cos(xhat(3)) 0;
        0 0 0 0 0 1;
        -inv(M)*[0 0 0 Xu+2*Xuu*abs(xhat(4)) -m22*xhat(6) -m22*xhat(5)-(m23+m32)*xhat(6);
        0 0 0 Yr*xhat(6)+m11*xhat(6) Yv+2*Yvv*abs(xhat(5)) Yr+m11*xhat(4); 
        0 0 0 m22*xhat(5)+((m23+m32)/2)*xhat(6)-m11*xhat(5) m22*xhat(4)+Nv-m11*xhat(4) ((m23+m32)/2)*xhat(4)+Nr+2*Nrr*abs(xhat(6))]]; % (replace with proper linearization of your dynamics if available)
    
    % State prediction and error covariance propagation
    Pmin  = FX*Pplus*FX' + QF*dt;
    Sigma = Pmin + RF*dt;   % measurement covariance
    
    % Kalman gain for state update
    KF = Pmin / Sigma;
    
    % Update state estimate
    xhat = xhat + KF*(y - xhat);
    Pplus = (eye(n_x) - KF)*Pmin;
    
    % Fault estimation using adaptive law:
    % Update fault estimation gain
    Upsilon = (eye(n_x) - KF)*FX*UpsilonPlus + (eye(n_x) - KF)*Psi;
    Omega   = FX*UpsilonPlus + Psi;
    Lambda  = inv(lambda*Sigma + Omega*S*Omega');
    Gamma   = S*Omega'*Lambda;
    S       = (1/lambda)*S - (1/lambda)*S*Omega'*Lambda*Omega*S;
    UpsilonPlus = Upsilon;
    
    % Update fault parameter estimate
    thetahat = thetahat + Gamma*(y - xhat);
    
    % Store values for plotting
    xArray(:,k)      = x;
    xhatArray(:,k)   = xhat;
    thetaArray(:,k)  = theta;
    thetahatArray(:,k) = thetahat;
    uArray(:,k)      = tau_c;
end

%% Plotting results
figure(1)
plot(xArray(1,:), xArray(2,:), 'k', 'LineWidth', 2); hold on;
plot(xhatArray(1,:), xhatArray(2,:), 'r--', 'LineWidth', 2);
legend('True Trajectory', 'AEKF Estimate');
xlabel('x (m)'); ylabel('y (m)'); grid on;
title('Otter Trajectory');

figure(2)
subplot(3,1,1)
plot(t, thetaArray(1,:), 'k', t, thetahatArray(1,:), 'r--', 'LineWidth', 2);
ylabel('\theta_1'); legend('True','Estimated'); grid on;
subplot(3,1,2)
plot(t, thetaArray(2,:), 'k', t, thetahatArray(2,:), 'r--', 'LineWidth', 2);
ylabel('\theta_2'); legend('True','Estimated'); grid on;
subplot(3,1,3)
plot(t, thetaArray(3,:), 'k', t, thetahatArray(3,:), 'r--', 'LineWidth', 2);
ylabel('\theta_3'); xlabel('Time (s)'); legend('True','Estimated'); grid on;
sgtitle('Actuator Fault Diagnosis');

figure(3)
plot(t, uArray(1,:), 'b', 'LineWidth', 2); hold on;
plot(t, uArray(2,:), 'r:', 'LineWidth', 2);
plot(t, uArray(3,:), 'g-.', 'LineWidth', 2);
legend('Control u_1','Control u_2','Control u_3');
xlabel('Time (s)'); ylabel('Control Input'); grid on;
title('Control Inputs');
