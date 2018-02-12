clear;

my_pi = 3.14159265358;

% -------- LQR design ---------

Q = [  10    0    0   ;
       0    10    0   ;
       0    0    5   ];
R = [  0.1    0   ;
       0    10   ];

% A = zeros(3);
% sq2 = sqrt(2);
% B = [sq2 0; sq2 0; 0 1];
% 
% [K,S,e] = lqr(A,B,Q,R);
% 
% K
% S
% e

x = sym('x', [3 1]);
u = sym('u', [2 1]);

f_original = [ u(1) * cos(x(3)) ; u(1) * sin(x(3)) ; u(2)];
eq_x = [12.5; 12.5; my_pi/4.0 ];

% f_original = [ - u(1) * cos(x(2)) ; -u(1) * sin(x(2)) / x(1) + u(2) ; u(1) * sin(x(2)) / x(1) ];
% eq_x = [0; 0; 0];

eq_u = [ 0.0; 0.0 ];

[A_l, B_l] = sym_linearization(f_original, x, u, eq_x, eq_u);

A_l
B_l

[K, S, ~] = lqr(A_l, B_l, Q, R);

K
S


% --------- run simulator -----------

robotCurrentPose = [[10.0 10.0] my_pi/4.0 + 0.3 ];

robotRadius = 0.4;
robot = ExampleHelperRobotSimulator('emptyMap',2);
robot.enableLaser(false);
robot.setRobotSize(robotRadius);
robot.showTrajectory(true);
robot.setRobotPose(robotCurrentPose);

controlRate = robotics.Rate(5);
robotGoal = [12.5 12.5];
angleGoal = my_pi/3.0;

distanceToGoal = norm(robotCurrentPose(1:2) - robotGoal);

theta = 0.0;

while (distanceToGoal > 0.1 || delta_theta > 0.1 )
    
    delta_x = robotCurrentPose(1) - eq_x(1);
    delta_y = robotCurrentPose(2) - eq_x(2);
    delta_theta = robotCurrentPose(3) - eq_x(3);
    
    % delta_x
    % delta_y
    
    U = -K * [delta_x; delta_y; delta_theta];
    
    % U
    
    theta = theta + U(2) / 2.0 * 0.12;
    
    drive(robot, U(1), theta);
    
    % disp(robotCurrentPose);

    waitfor(controlRate);
    
    % Re-compute the distance to the goal
    robotCurrentPose = robot.getRobotPose;
    distanceToGoal = norm(robotCurrentPose(1:2) - robotGoal);
    delta_theta = abs(robotCurrentPose(3) - my_pi/4.0);
    
    % distanceToGoal
    % delta_theta
    
end
drive(robot, 0, 0);
robotCurrentPose = robot.getRobotPose;
% distanceToGoal
% delta_theta