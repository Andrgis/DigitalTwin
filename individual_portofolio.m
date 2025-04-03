%% Individual Portofolio - Otter
clc; clear; clf;

%%  Formulate the discrete-time model of the Otte
% Parameters from Table 1
Xu = -0.7225; Xuu = -1.3274; m = 23.8; X_udot = -2.0;
Yv = -0.8612; Yvv = -36.2823; Yv_dot = -10.0; xg = 0.046;
Yr = 0.1079; Nrr = -1; Yr_dot = 0.0; Nv_dot = 0.0;
Nv = 0.1052; Nr = -0.5; Izz = 1.76; Nr_dot = -1.0;


% Variables:
dt = 0.0005;
T = 15;
t = dt:dt:T;

% Noise
QF = 10*eye(6);
RF = 40*eye(6);

%Initial:
    x = 0;     y = 0;   yaw = 0;
    u = 2;     v = 1;     r = 0;
tau_u = -30; tau_v = 0; tau_r = 0;

eta   = [x, y, yaw]';
nu    = [u, v, r]';
tau_c = [tau_u, tau_v, tau_r]';
tau_nf = tau_c;
X     = [eta; nu];
X_nf  = X;
xhat = [2,2,0,1,1,0]';
theta = [0,0,0]';
thetahat = [0,0,0]';


%% Tuning parameters
lambda      = 0.995;%0.995;
a           = 0.999;%0.999;

%% Mass matrix

A = eye(6);
% Mass matrix, Coriolis, and damping (as in your simulation code)
M = [ m - X_udot,           0,            0;
          0,     m - Yv_dot,  m*xg - Yr_dot;
          0,     m*xg - Nv_dot, Izz - Nr_dot ];
        
m11 = M(1,1);
m22 = M(2,2);
m23 = M(2,3);
m32 = M(3,2);
            
B_c = [zeros(3); 
       inv(M)*eye(3)];   
   
Psi_d = -dt*B_c*diag(tau_c);
%% Estimation paramaters
Pplus  = eye(rank(A));
S      = 0.1 * eye(3);  % Use a 3x3 gain matching fault dimension
UpsilonPlus = 0*B_c;

C = eye(6);

x_nf_vec  = zeros(length(X_nf) ,length(t));
eta_vec   = zeros(length(eta)  ,length(t));
nu_vec    = zeros(length(nu)   ,length(t));
tau_c_vec = zeros(length(tau_c),length(t));
theta_vec = zeros(length(theta),length(t));
thetahat_vec = zeros(length(thetahat),length(t));
xhat_vec  = zeros(length(xhat) ,length(t));

%% Simulation no fault
for i = 1:length(t)
    if i==round(0.25*length(t))
        tau_nf = [200 50 0.2]';
    end
    
    if i==round(0.5*length(t))
        tau_nf = [50 -40 -0.5]';        
    end
    
    % Rotation matrix
    R = [cos(X_nf(3)), -sin(X_nf(3)), 0; %Correct
         sin(X_nf(3)),  cos(X_nf(3)), 0;
                0,         0, 1];

    % Corolis Matrix
    C_M = [0, 0,        -M(2,2)*X_nf(5)-0.5*(M(2,3)+M(3,2))*X_nf(6); %Correct
           0, 0,                               M(1,1)*X_nf(4);
           M(2,2)*X_nf(5)+0.5*(M(2,3)+M(3,2))*X_nf(6), -M(1,1)*X_nf(4), 0];

    % Damping matrix
    D_M = -[Xu+Xuu*abs(X_nf(4)),             0,            0; %Correct
                        0, Yv+Yvv*abs(X_nf(5)),           Yr;
                        0,            Nv, Nr+Nrr*abs(X_nf(6))];
                    
    A_c = [zeros(3),      R;
           zeros(3), -inv(M)*(C_M+D_M)];
    
    % Discreetizing
    A_d = eye(size(A_c)) + dt*A_c;
    B_d = dt*B_c;
    
    X_nf = A_d*X_nf + B_d*tau_nf + QF*dt*randn(size(X_nf));
    
    x_nf_vec(:,i)  = X_nf;
end

%% Simulation fault
for i = 1:length(t)

    % Rotation matrix
    R = [cos(eta(3)), -sin(eta(3)), 0; %Correct
         sin(eta(3)),  cos(eta(3)), 0;
                0,         0, 1];

    % Corolis Matrix
    C_M = [0, 0,        -M(2,2)*nu(2)-0.5*(M(2,3)+M(3,2))*nu(3); %Correct
           0, 0,                               M(1,1)*nu(1);
           M(2,2)*nu(2)+0.5*(M(2,3)+M(3,2))*nu(3), -M(1,1)*nu(1), 0];

    % Damping matrix
    D_M = -[Xu+Xuu*abs(nu(1)),             0,            0; %Correct
                        0, Yv+Yvv*abs(nu(2)),           Yr;
                        0,            Nv, Nr+Nrr*abs(nu(3))];
                    
    A_c = [zeros(3),      R;
           zeros(3), -inv(M)*(C_M+D_M)];
    
    % Discreetizing
    A_d = eye(size(A_c)) + dt*A_c;
    B_d = dt*B_c;
    
    if i==round(0.25*length(t))
        tau_c = [200 50 0.2]';
    end
    if i==round(0.5*length(t))
        theta = [0 0.041 0]';
        tau_c = [50 -40 -0.5]';        
    end
    if i==round(0.75*length(t))
        theta = [0.087 0.041 0.2]';       
    end
    
    % Fault design
    Psi_d = -B_d*diag(tau_c);

    % Updating variables
%     eta = eta + dt * R*nu;
%     nu = nu + dt * inv(M)*(-(C_M + D_M)*nu + eye(3)*tau_c);
   %X = State development + Control influence 
   %  + Actuator fault + model uncertainty
    X = A_d*X + B_d*tau_c + Psi_d*theta + QF*dt*randn(size(X));
    eta = X(1:3);
    nu = X(4:6);
    y = C*X + RF*dt*randn(size(X));
    
    % Adaptive Extended Kalman filter
    %Jacobian 
    sig2 = 0.5 * (m23 + m32);
    sig1 = xhat(5) * m22 + xhat(6) * sig2;

    nudot = -inv(M)*[Xu+2*Xuu*abs(xhat(4)), xhat(6)*m22, sig1+xhat(6)*sig2;
                 -xhat(6)*m11, Yv+2*Yvv*abs(xhat(5)), Yr-xhat(4)*m11;
                  xhat(5)*m11-sig1, Nv+xhat(4)*m11-xhat(4)*m22, Nr+2*Nrr*abs(xhat(6))-xhat(4)*sig2];

    FX = A + dt * [0, 0, -xhat(5)*cos(xhat(3))-xhat(4)*sin(xhat(3)), cos(xhat(3)), -sin(xhat(3)), 0;
                   0, 0,  xhat(4)*cos(xhat(3))-xhat(5)*sin(xhat(3)), sin(xhat(3)),  cos(xhat(3)), 0;
                   0, 0,  0, 0, 0, 1;
                   zeros(3), nudot];

    Pmin  = FX*Pplus*FX'+QF;
    Sigma = C*Pmin*C'+RF;
    KF    = Pmin*C'*inv(Sigma);
    Pplus = (eye(rank(A))-KF*C)*Pmin;
    
    ytilde = y-C*xhat;
    QF    = a*QF + (1-a)*(KF*(ytilde*ytilde')*KF');    
    RF    = a*RF + (1-a)*(ytilde*ytilde'+C*Pmin*C');

    Upsilon = (eye(rank(A))-KF*C)*FX*UpsilonPlus+(eye(rank(A))-KF*C)*Psi_d;
    Omega   = C*FX*UpsilonPlus+C*Psi_d;
    Lambda  = inv(lambda*Sigma+Omega*S*Omega');
    Gamma   = S*Omega'*Lambda;
    S       = (1/lambda)*S-(1/lambda)*S*Omega'*Lambda*Omega*S;
    UpsilonPlus = Upsilon;
    
    thetahat  = thetahat + Gamma*(y-C*xhat);
    xhat      = A_d*xhat + B_d*tau_c + Psi_d*thetahat + QF*dt*randn(size(X)) +KF*(y-C*xhat)+Upsilon*Gamma*(y-C*xhat);
    
    % Storing values
    eta_vec(:,i)   = X(1:3);
    nu_vec(:,i)    = X(4:6);
    tau_c_vec(:,i) = tau_c;
    xhat_vec(:,i)  = xhat;
    theta_vec(:,i) = theta;
    thetahat_vec(:,i) = clamp(thetahat);
end


%% Plotting
% Plot 1 - position simulation and approximation
figure(1)
clf;
hold on
plot(eta_vec(1,:), eta_vec(2,:), 'k', 'LineWidth', 3)
plot(xhat_vec(1,:),xhat_vec(2,:), 'r:', 'LineWidth', 3)
% plot(x_nf_vec(1,:),x_nf_vec(2,:), 'g', 'LineWidth', 1.5)
plot(eta_vec(1,round(0.25*T/dt)), eta_vec(2,round(0.25*T/dt)), "gpentagram", "LineWidth", 3)
plot(eta_vec(1,round(0.5*T/dt)), eta_vec(2,round(0.5*T/dt)), "bpentagram", "LineWidth", 3)
plot(eta_vec(1,round(0.75*T/dt)), eta_vec(2,round(0.75*T/dt)), "rpentagram", "LineWidth", 3)
grid minor
ylabel('y (m)','FontSize',12)
xlabel('x (m)','FontSize',12)
% legend("simulation", "approximation", "No fault")
legend("simulation", "approximation", "\tau 1", "\tau 2, fault 1", "fault 2")

%% Plot 2 - Vessel frame parameters


figure(2);
clf;
subplot(3,1,1)
plot(t,nu_vec(1,:), 'k', 'LineWidth', 3)
hold on
plot(t,xhat_vec(4,:), 'r:', 'LineWidth', 3)
plot(0.25*T, nu_vec(1,round(0.25*T/dt)), "gpentagram", "LineWidth", 3)
plot(0.5*T, nu_vec(1,round(0.5*T/dt)), "bpentagram", "LineWidth", 3)
plot(0.75*T, nu_vec(1,round(0.75*T/dt)), "rpentagram", "LineWidth", 3)
grid on;
grid minor
ylabel('u (m/s)','FontSize',12)
set(gca,'FontSize',12)
subplot(3,1,2)
plot(t,nu_vec(2,:), 'k', 'LineWidth', 3)
hold on
plot(t,xhat_vec(5,:), 'r:', 'LineWidth', 3)
plot(0.25*T, nu_vec(2,round(0.25*T/dt)), "gpentagram", "LineWidth", 3)
plot(0.5*T, nu_vec(2,round(0.5*T/dt)), "bpentagram", "LineWidth", 3)
plot(0.75*T, nu_vec(2,round(0.75*T/dt)), "rpentagram", "LineWidth", 3)
grid on;
grid minor
ylabel('v (m/s)','FontSize',12)
subplot(3,1,3)
plot(t,nu_vec(3,:), 'k', 'LineWidth', 3)
hold on
plot(t,xhat_vec(6,:), 'r:', 'LineWidth', 3)
plot(0.25*T, nu_vec(3,round(0.25*T/dt)), "gpentagram", "LineWidth", 3)
plot(0.5*T, nu_vec(3,round(0.5*T/dt)), "bpentagram", "LineWidth", 3)
plot(0.75*T, nu_vec(3,round(0.75*T/dt)), "rpentagram", "LineWidth", 3)
grid on;
grid minor
ylabel('r (rads/s)','FontSize',12)
set(gca,'FontSize',12)
xlabel('Time (s)','FontSize',12)
set(gca,'FontSize',12)
h1 = legend('True State','AEKF', "\tau 1", "\tau 2, fault 1", "fault 2",'FontSize',12);
set(h1, 'Position', [0.7, 0.8, .1, .1])

%% Plot - Fault parameter estimation

figure(3)
clf;
subplot(3,1,1)
plot(t,theta_vec(1,:), 'k', 'LineWidth', 3)
hold on;
plot(t,thetahat_vec(1,:), 'r:', 'LineWidth', 3)
plot(0.25*T, theta_vec(1,round(0.25*T/dt)), "gpentagram", "LineWidth", 3)
plot(0.5*T, theta_vec(1,round(0.5*T/dt)), "bpentagram", "LineWidth", 3)
plot(0.75*T, theta_vec(1,round(0.75*T/dt)), "rpentagram", "LineWidth", 3)
ylabel('\theta_1','FontSize',12)
xlabel('Time (s)','FontSize',12)
grid on
grid minor
set(gca,'FontSize',12)
subplot(3,1,2)
plot(t,theta_vec(2,:), 'k', 'LineWidth', 3)
hold on;
plot(t,thetahat_vec(2,:), 'r:', 'LineWidth', 3)
plot(0.25*T, theta_vec(2,round(0.25*T/dt)), "gpentagram", "LineWidth", 3)
plot(0.5*T, theta_vec(2,round(0.5*T/dt)), "bpentagram", "LineWidth", 3)
plot(0.75*T, theta_vec(2,round(0.75*T/dt)), "rpentagram", "LineWidth", 3)
ylabel('\theta_2','FontSize',12)
xlabel('Time (s)','FontSize',12)
grid on
grid minor
set(gca,'FontSize',12)
subplot(3,1,3)
plot(t,theta_vec(3,:), 'k', 'LineWidth', 3)
hold on;
plot(t,thetahat_vec(3,:), 'r:', 'LineWidth', 3)
plot(0.25*T, theta_vec(3,round(0.25*T/dt)), "gpentagram", "LineWidth", 3)
plot(0.5*T, theta_vec(3,round(0.5*T/dt)), "bpentagram", "LineWidth", 3)
plot(0.75*T, theta_vec(3,round(0.75*T/dt)), "rpentagram", "LineWidth", 3)
legend('True \theta','Estimated \theta', "\tau 1", "\tau 2, fault 1", "fault 2",'FontSize',12);
ylabel('\theta_3','FontSize',12)
xlabel('Time (s)','FontSize',12)
grid on
grid minor
set(gca,'FontSize',12)

%% Plot - controll input

figure(4)
clf;
subplot(3,1,1)
plot(t,tau_c_vec(1,:), 'k', 'LineWidth', 3)
hold on;
ylabel('\tau_u','FontSize',12)
xlabel('Time (s)','FontSize',12)
grid on
grid minor
set(gca,'FontSize',12)
subplot(3,1,2)
plot(t,tau_c_vec(2,:), 'k', 'LineWidth', 3)
hold on;
ylabel('\tau_v','FontSize',12)
xlabel('Time (s)','FontSize',12)
grid on
grid minor
set(gca,'FontSize',12)
subplot(3,1,3)
plot(t,tau_c_vec(3,:), 'k', 'LineWidth', 3)
hold on;
ylabel('\tau_r','FontSize',12)
% legend('True \theta','Estimated \theta', "\tau 1", "\tau 2, fault 1", "fault 2",'FontSize',12);
xlabel('Time (s)','FontSize',12)
grid on
grid minor
set(gca,'FontSize',12)

%% functions

function clamped_arr = clamp(arr)
    clamped_arr = arr;
    for i=1:length(arr)
        x = arr(i);
        if x>1
            x = 1;
        end
        if x<-1
            x = -1;
        end
        clamped_arr(i) = x;
    end
end