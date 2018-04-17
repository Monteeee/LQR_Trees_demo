%% basic state model and time invariant LQR around the final state

clear;
mypi = 3.14159265358;
% tolerance (slack variable to transfer negativity to non-positivity)
epsilon = 1e-7;
% initial value of rho
rho = 10.0;

verbose = 4;

Q = [  10   0   ;
       0    1   ];
R = 15;

x = sym('x', [2 1]);
u = sym('u', [1 1]);
Nx = length(x);
x_bar = sym('x_', [Nx, 1]);
x_ = x_bar;

m = 1.0;
g = 9.8;
l = 0.5;
I = m*l^2;
b = 0.1;

f_original = [ x(2) ; 
               -1/I * ( b * x(2) + m*g*l*sin(x(1)) ) + u(1)/I ];

eq_x = [mypi; 0.0];
eq_u = [ 0.0 ];

[A_l, B_l] = sym_linearization(f_original, x, u, eq_x, eq_u);

A_l
B_l

N = zeros(2 , 1);

[K, S, ~] = lqr(A_l, B_l, Q, R);

K
S


%% simulation of TI LQR control

if verbose == 1 || verbose == 2

    x_init = [mypi - 0.2; 0.0];

    controlRate = robotics.Rate(10);

    F = 0;
    dt = 1.0/10;

    A = A_l;
    B = B_l;

    iter = 50;
    recordl = zeros(iter, 1);
    recordv = zeros(iter, 1);
    recordu = zeros(iter, 1);

    delta = x_init - eq_x;
    for i=1:iter

        F = - K * delta;

        if F > 3
            F = 3;
        elseif F < -3
            F = -3;
        end

        d1 = delta(2);
        d2 = -1/I * ( b * delta(2) + m*g*l*sin(delta(1)) ) + F/I;
        dl = [d1;d2];

        recordl(i) = delta(1);
        recordv(i) = delta(2);
        recordu(i) = F;

        waitfor(controlRate);
        delta = delta + dl.*dt;
    end

    if verbose == 2
        plot(recordl,'DisplayName','recordl');
        hold on;
        plot(recordu,'DisplayName','recordu');
        plot(recordv,'DisplayName','recordv');
        hold off;
    end
end
%% verification of the basin of attraction of TI LQR

if verbose == 3 || verbose == 4

    % define f_hat

    % closed loop function
    f_cl = toCL(f_original, x, u, eq_x, eq_u, K);

    f_cl

    % taylor expansion
    f_cl_t = taylor(f_cl, x.', eq_x.', 'order', 4);

    temp = vpa(f_cl_t, 5)

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

    f_hat

    % define dJ*_hat(x_)
    dJ = 2 .* (x_.') * S * f_hat;

    % define J*(x_)
    J = (x_.') * S * x_;

    % define norm part
    nor22 = epsilon * ( x_(1)^2 + x_(2)^2 );

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

%% plot out current region of attraction of this controller


if verbose == 4

    ttt = vpa(simplify(J), 5);

    x1gv = linspace(-mypi, mypi, 1000); % grid vector for theta
    x2gv = linspace(-4*mypi, 4*mypi, 1000); % grid vector for d(theta)/dt

    [vx, vy] = meshgrid(x1gv, x2gv);

    % tt = subs(subs(J__, x_(1), vx), x_(3), vy);

    % f_(x_) = ttt;
    % tt = double(f_(x1gv, x2gv));
    
    tt = 174.14*vx.^2 + 74.007*vx.*vy + 8.019*vy.^2;

    condition1 = tt <= rho;
    condition2 = tt >= 0;
    output = ones(length(x1gv), length(x2gv)); % Initialize to 1
    output(~(condition1&condition2)) = 0; % Zero out coordinates not meeting conditions.
    imshow(output, 'xdata', x1gv, 'ydata', x2gv); % Display
    axis on;
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ------   --------
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
end

