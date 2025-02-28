%% Research code by Agus Hasan

clear;
clc;

%% Time horizon
tf  = 20;
dt  = 0.001;
t   = dt:dt:tf;

%% System description
A = eye(4);
B = dt*[0 0 1 0;0 0 0 1]';
C = eye(4);

%% Noise
QF = 0.01*eye(rank(A));
RF = 0.04*eye(rank(C));

%% Initialization
x        = [0;0;0;0];

%% Control
u           = [1 1]';

%% For plotting
uArray     = [];
xArray     = [];

%% Simulation
for i=1:(tf/dt)
    
    if i>5000
        u = [2 0.5]';
    end
    
    if i>10000
        u = -[-0.5 0.4]';        
    end
    
   uArray         = [uArray u];
   xArray         = [xArray x];
      
   x = A*x+dt*[x(3)*sin(x(4)) x(3)*cos(x(4)) 0 0]'+B*u+QF*dt*randn(4,1);
   y = C*x+RF*dt*randn(4,1);
end

figure(1);
plot(xArray(1,:),xArray(2,:), 'k', 'LineWidth', 3)
grid on;
grid minor
ylabel('y (m)','FontSize',12)
xlabel('x (m)','FontSize',12)
set(gca,'FontSize',12)
legend('State','FontSize',12);

figure(2);
subplot(2,1,1)
plot(t,xArray(3,:), 'k', 'LineWidth', 3)
grid on;
grid minor
ylabel('v (m/s)','FontSize',12)
set(gca,'FontSize',12)
subplot(2,1,2)
plot(t,xArray(4,:), 'k', 'LineWidth', 3)
grid on;
grid minor
ylabel('\psi (rad)','FontSize',12)
xlabel('Time (s)','FontSize',12)
set(gca,'FontSize',12)
h1 = legend('State','FontSize',12);
set(h1, 'Position', [0.7, 0.8, .1, .1])

figure(3);
plot(t,uArray(1,:), 'b', 'LineWidth', 3)
hold on;
plot(t,uArray(2,:), ':r', 'LineWidth', 3)
grid on;
grid minor
ylabel('u','FontSize',12)
xlabel('Time (s)','FontSize',12)
legend('a','r')
set(gca,'FontSize',12)
































