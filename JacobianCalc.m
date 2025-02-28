% Parameters from Table 1
Xu = -0.7225; Xuu = -1.3274; m = 23.8; X_udot = -2.0;
Yv = -0.8612; Yvv = -36.2823; Yv_dot = -10.0; xg = 0.046;
Yr = 0.1079; Nrr = -1; Yr_dot = 0.0; Nv_dot = 0.0;
Nv = 0.1052; Nr = -0.5; Izz = 1.76; Nr_dot = -1.0;

M = [m-X_udot,           0,           0; %Correct
                0,    m-Yv_dot, m*xg-Yr_dot;
                0, m*xg-Nv_dot,  Izz-Nr_dot];
            
syms x y t u v r
eta = [x, y, t]';
nu = [u, v, r]';
X = [eta;nu];

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
       zeros(3), -(C_M+D_M)/M];
   
J = jacobian(A_c*X,X')