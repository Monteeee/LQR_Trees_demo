%% basic state model and time invariant LQR around the final state
% ------------------------------------------------------------------------

clear;
mypi = 3.14159265358;
% tolerance (slack variable to transfer negativity to non-positivity)
epsilon = 1e-10;
% initial value of rho
rho = 0.43;

verbose = 4;

Q = [  5    0    0    0;
       0    5    0    0;
       0    0    10   0;
       0    0    0    10];
R = 10;

x = sym('x', [4 1]);
u = sym('u', [1 1]);
Nx = length(x);
Nu = length(u);
x_bar = sym('x_', [Nx, 1]);
x_ = x_bar;

m1 = 2.0;
m2 = 0.5;
g = 9.81;
l = 1.0;
d = 1.0;

f_original = [ x(2) ;
               ( l*m2*sin(x(3))*x(4)^2 + u(1) + m2*g*cos(x(3))*sin(x(3)) ) / ( m1 + m2*(1 - (cos(x(3)))^2) );
               x(4) ;
               - ( l*m2*cos(x(3))*sin(x(3))*x(4)^2 + u(1)*cos(x(3)) + (m1+m2)*g*sin(x(3)) ) / ( l*m1 + l*m2*(1 - (cos(x(3)))^2) ) ];

f_dynamics = @(x, u) [ x(2) ;
               ( l*m2*sin(x(3))*x(4)^2 + u(1) + m2*g*cos(x(3))*sin(x(3)) ) / ( m1 + m2*(1 - (cos(x(3)))^2) );
               x(4) ;
               - ( l*m2*cos(x(3))*sin(x(3))*x(4)^2 + u(1)*cos(x(3)) + (m1+m2)*g*sin(x(3)) ) / ( l*m1 + l*m2*(1 - (cos(x(3)))^2) ) ];
           
eq_x = [d ; 0.0; mypi; 0.0];
eq_u = [ 0.0 ];

init_x = [d - 0.01; 0.0; mypi - 0.01; 0.0];
init_u = [0.0];

[A_l, B_l] = sym_linearization(f_original, x, u, eq_x, eq_u);

A_l
B_l

[K, S, ~] = lqr(A_l, B_l, Q, R);

K
S


%% simulation of TI LQR control'
% --------------------------------------

if verbose == 1 || verbose == 2

    controlRate = robotics.Rate(1000);

    F = 0;
    dt = 1.0/1000;

    iter = 100;
    record1 = zeros(iter, 1);
    record2 = zeros(iter, 1);
    recordu = zeros(iter, 1);

    delta = init_x - eq_x;
    
    for i=1:iter

        F = - K * delta;

%         if F > 3
%             F = 3;
%         elseif F < -3
%             F = -3;
%         end

        dl = [delta(2);
              ( l*m2*sin(delta(3))*delta(4)^2 + F + m2*g*cos(delta(3))*sin(delta(3)) ) / ( m1 + m2*(1 - (cos(delta(3)))^2) );
              delta(4);
               -( l*m2*cos(delta(3))*sin(delta(3))*delta(4)^2 + F*cos(delta(3)) + (m1+m2)*g*sin(delta(3)) ) / ( l*m1 + l*m2*(1 - (cos(delta(3)))^2) )
              ];
        
%         dl = A_l * delta + B_l * F;

        record1(i) = delta(1);
        record2(i) = delta(3);
        recordu(i) = F;

        waitfor(controlRate);
        delta = delta + dl.*dt;
    end

    if verbose == 2
        plot(record1,'DisplayName','recordl');
        hold on;
        plot(recordu,'DisplayName','recordu');
        plot(record2,'DisplayName','recordv');
        hold off;
    end
end

%% verification of the basin of attraction of TI LQR
%  ----------------------------------------------------------

if verbose == 3 || verbose == 4

    % define f_hat

    % closed loop function
    f_cl = toCL(f_original, x, u, eq_x, eq_u, K);

    f_cl;

    % taylor expansion
    f_cl_t = taylor(f_cl, x.', eq_x.', 'order', 4);

    temp = vpa(f_cl_t, 5);

    % substitute the state x with error state x_
    for k = 1:length(x)
        f_cl_t = subs(f_cl_t, x(k), x_bar(k)+eq_x(k));
    end

    for i=1:length(f_cl_t)
        [C, T] = coeffs(f_cl_t(i));
        C(abs(C) < 1e-7) = 0;
        f_cl_t(i) = dot(C, T);
    end

    f_hat = vpa(f_cl_t, 5);

    f_hat;

    % define dJ*_hat(x_)
    dJ = 2 .* (x_.') * S * f_hat;
    dJ__ = vpa(dJ, 5);
    
    % define J*(x_)
    J = (x_.') * S * x_;
    J__ = vpa(J, 5);
    
    % define norm part
    nor22 = epsilon * ( x_(1)^2 + x_(2)^2 + x_(3)^2 + x_(4)^2 );

    % define a sos program
    Program1 = sosprogram(x_);

    % define h(x_) as sums of squares
    [Program1, h] = sossosvar(Program1, x_);

    % add inequality constraint
    Program1 = sosineq(Program1, -dJ - h*(rho - J) - nor22);

    % set solver option

    Program1 = sossolve(Program1);

    SOLV = sosgetsol(Program1, h);

    disp(SOLV)
end


%% try drawing the elliptical region of attraction
% -------------------------------------------------

% to make a 2-D and illustrative plot, assume dx/dt and dtheta/dt are
% both 0 (equilibrium state) in initial condition.

if verbose == 4
    J__ = subs(J__, x_(2), 0.0);
    J__ = subs(J__, x_(4), 0.0);

    ttt = vpa(simplify(J__), 5);

    x1gv = linspace(-0.3, 0.3, 1000); % grid vector for position
    x3gv = linspace(-mypi/8.0, mypi/8.0, 1000); 

    [vx, vy] = meshgrid(x1gv, x3gv);

    % tt = subs(subs(J__, x_(1), vx), x_(3), vy);

    tt = 13.245*vx.^2 - 139.1*vx.*vy + 1547.6*vy.^2;

    condition1 = tt <= rho;
    condition2 = tt >= 0;
    output = ones(length(x1gv), length(x3gv)); % Initialize to 1
    output(~(condition1&condition2)) = 0; % Zero out coordinates not meeting conditions.
    imshow(output, 'xdata', x1gv, 'ydata', x3gv); % Display
    axis on;
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ------  the estimated controllable region is extremely small --------
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
end


%% some trajectory planning example

if verbose == 5

    % in real case this will be the sampled point
    startPoint = [d - 0.2; 0.0; mypi - 0.1; 0.0]; % [x, dx/dt, theta, dtheta/dt]

    % assum we have already got some point in the tree
    tree_point = [d,        0.0, mypi,        0.0;
                  d - 0.05, 0.0, mypi - 0.05, 0.0;
                  d - 0.1, 0.08, mypi + 0.03, -0.01 ]';

    finishPoint = tree_point(:, 1);

    finishPoint_upp = max(tree_point, [], 2);
    finishPoint_low = min(tree_point, [], 2);

    max_force = 50.0;

    % add in dynamics, objective and finalconstraint (tree point constraint)
    problem.func.dynamics = @(t,x,u)( dynamics(x,u) );
    problem.func.pathObj = @(t,x,u)( objective(x,u) );
    problem.func.bndCst = @(t0,x0,tF,xF)( finalconstraint(x0,xF) );

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %                 Set up bounds on state and control                      %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    problem.bounds.initialTime.low = 0;
    problem.bounds.initialTime.upp = 0;
    problem.bounds.finalTime.low = 5;
    problem.bounds.finalTime.upp = 10;

    problem.bounds.state.low = [-5; -10; -2*mypi; -2*mypi ]; % bound along trajectory
    problem.bounds.state.upp = [ 5;  10;  2*mypi;  2*mypi ];

    problem.bounds.initialState.low = startPoint;
    problem.bounds.initialState.upp = startPoint;

    problem.bounds.finalState.low = finishPoint_low;
    problem.bounds.finalState.upp = finishPoint_upp;

    problem.bounds.control.low = -max_force;
    problem.bounds.control.upp = max_force;

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %                 Initialize trajectory with guess                        %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % Car travels at a speed of one, and drives in a straight line from start
    % to finish point.

    del = finishPoint - startPoint;  % vector from start to finish

    problem.guess.time = [0, norm(del)];   % time = distance/speed
    problem.guess.state = [ startPoint, finishPoint];
    problem.guess.control = [0, 0];  

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %                      Options for Transcription                          %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    problem.options.nlpOpt = optimset(...
        'display','iter',...
        'MaxFunEval',1e5,...
        'tolFun',1e-6);

    problem.options.method = 'hermiteSimpson';
    problem.options.hermiteSimpson.nSegment = 25;

    % problem.options.method = 'gpops';

    % problem.options.method = 'trapezoid';
    % problem.options.trapezoid.nGrid = 15;

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %                            Solve!                                       %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

    soln = optimTraj(problem);

    t = linspace(soln.grid.time(1), soln.grid.time(end), 150);
    z_data = soln.interp.state(t);
    x_data = z_data(1,:);
    y_data = z_data(3,:);
    u_data = soln.interp.control(t);
    
    plot(x_data,'DisplayName','x_data');hold on;plot(y_data,'DisplayName','y_data');hold off;
    
    
end


%% to solve the matrix ODE for S(t)

x_t = sym('x0_', [Nx, 1]);
u_t = sym('u0_', [Nu, 1]);

[A_t, B_t] = sym_linearization(f_original, x, u, x_t, u_t);

% save the data so that in Riccarti function I can use them without reading
% them in every function call
% save('opti_traj.mat', 'soln', 'A_t', 'B_t', 'x_t', 'u_t', 'Q', 'R');

disp("done!");


%% the end
disp("done!");