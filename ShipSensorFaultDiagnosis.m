%% Research code by Agus Hasan

clear;
clc;

%% Time horizon
tf  = 40;
dt  = 0.01;
t   = dt:dt:tf;

%% System description
A = eye(4);
B = [0 0;0 0;dt 0;0 dt];
C = [0 0 1 0; 0 0 0 1];

%% Filtering matrix
Af  = 0.5*eye(2);

%% Augmented system
Ab  = [A 0*eye(4,2);
       Af*dt*C eye(2)-Af*dt];
Bb  = [B;0*eye(2)];
Cb  = [0 0 0 0 1 0;
      0 0 0 0 0 1];

%% Noise
QF          = 0.001*eye(rank(Ab));
RF          = 0.001*eye(2);

%% Control
u = [1 0.1]';

%% Fault term
Phi = [0*eye(4,2);
       Af*dt];

%% Initialization
x             = [0 0 0 0 0 0]';
xhat          = [0 0 0 0 0 0]';
theta         = [0.0;0.0];
thetahat      = [0;0];

%% Estimation parameters
Pplus       = 1*eye(rank(Ab));
S           = 0.65*eye(2);
UpsilonPlus = 0*eye(6,2);

%% Tuning parameters
a           = 0.999;
lambda      = 0.995;

%% For plotting
xArray        = [];
xhatArray     = [];
thetaArray    = [];
thetahatArray = [];

%% Simulation

for k = 1:(tf/dt)
   
    if k>1000
        theta       = [0.98;0.31];
    end

    xArray          = [xArray x];
    xhatArray       = [xhatArray xhat];
    thetaArray      = [thetaArray theta];
    thetahatArray   = [thetahatArray thetahat];
    
    % Simulating the system
    x = Ab*x+dt*[x(3)*cos(x(4));x(3)*sin(x(4));0;0;0;0]+Bb*u+Phi*theta;
    % Taking measurement
    y = Cb*x;
    
    FX = Ab+dt*[0 0 cos(xhat(4)) -xhat(3)*sin(xhat(4)) 0 0; 0 0 sin(xhat(4)) xhat(3)*cos(xhat(4)) 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];

    % Estimation using observer
     Pmin  = FX*Pplus*FX'+QF;
     Sigma = Cb*Pmin*Cb'+RF;
     KF    = Pmin*Cb'*inv(Sigma);
     Pplus = (eye(rank(Ab))-KF*Cb)*Pmin;
     
     ytilde = y-Cb*xhat;
     QF    = a*QF + (1-a)*(KF*(ytilde*ytilde')*KF');    
     RF    = a*RF + (1-a)*(ytilde*ytilde'+Cb*Pmin*Cb');
 
    Upsilon = (eye(rank(Ab))-KF*Cb)*FX*UpsilonPlus+(eye(rank(Ab))-KF*Cb)*Phi;
    Omega   = Cb*FX*UpsilonPlus+Cb*Phi;
    Lambda  = inv(lambda*Sigma+Omega*S*Omega');
    Gamma   = S*Omega'*Lambda;
    S       = (1/lambda)*S-(1/lambda)*S*Omega'*Lambda*Omega*S;
    UpsilonPlus = Upsilon;
    
    thetahat  = thetahat + Gamma*(y-Cb*xhat);
    xhat      = Ab*xhat+dt*[xhat(3)*cos(xhat(4)) xhat(3)*sin(xhat(4)) 0 0 0 0]'+Bb*u+Phi*thetahat+KF*(y-Cb*xhat)+Upsilon*Gamma*(y-Cb*xhat);
end

%% Plotting
figure(1)
plot(xArray(1,:),xArray(2,:),'-','LineWidth',3)
hold on;
plot(xhatArray(1,:),xhatArray(2,:),':','LineWidth',3)
set(gca,'color','white','FontSize',14)
grid on;
grid minor;
ylabel('position (m)')
legend('actual','estimated')
xlabel('time (s)')

figure(2)
subplot(2,1,1)
plot(t,thetaArray(1,:),'-','LineWidth',3)
hold on;
plot(t,thetahatArray(1,:),':','LineWidth',3)
set(gca,'color','white','FontSize',14)
grid on;
grid minor;
ylabel('theta_1')
legend('actual','estimated')
subplot(2,1,2)
plot(t,thetaArray(2,:),'-','LineWidth',3)
hold on;
plot(t,thetahatArray(2,:),':','LineWidth',3)
set(gca,'color','white','FontSize',14)
grid on;
grid minor;
ylabel('theta_2')
xlabel('time (s)')