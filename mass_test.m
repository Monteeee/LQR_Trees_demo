
controlRate = robotics.Rate(10);

%% system dynamics


x = sym('x', [2 1]);
u = sym('u', [1 1]);
x_ = sym('x_', [2 1]);
u_ = sym('u_', [1 1]);

m = 1;
k = 0.5;

f_original = [ x(2) ; -k/m * x(1) + 1/m * u(1)];

init = [3;0];

l = init;

A = [0 1; -k/m 0];
B = [0 ; 1/m];

F = 0;
dt = 1.0/10;

recordl = zeros(100, 1);
recordv = zeros(100, 1);


%% LQR controller

eq_x = [2.0; 0];
eq_u = [eq_x(1)*k];

[A_lin, B_lin] = sym_linearization(f_original, x, u, eq_x, eq_u);

Q = [  5    0   ;
       0    5  ];
R = [  1  ];

[K, S, ~] = lqr(A_lin, B_lin, Q, R);

%% simulation

for i=1:100
    delta = l - eq_x;

    F = - K * delta;
    
    dl = A * l + B .* F;
    
    recordl(i) = l(1);
    recordv(i) = l(2);
    
    waitfor(controlRate);
    l = l + dl.*dt;
end

plot([1:100], recordl)

%% SOS verification

mypi = 3.14159265358;

% tolerance 
epsilon = 0.0001;
% initial value of rho
rho = 10000.0;

% define f_hat

% closed loop function
f_cl = toCL(f_original, x, u, eq_x, eq_u, K);

% taylor expansion
f_cl_t = symbolic_Taylor(f_cl, x, eq_x);

for i=1:length(f_cl_t)
    [C, T] = coeffs(f_cl_t(i));
    C(abs(C) < 1e-7) = 0;
    f_cl_t(i) = dot(C, T);
end
f_hat = vpa(f_cl_t, 5)

% define dJ*_hat(x_)
dJ = 2 .* (x_.') * S * f_hat;

% define J*(x_)
J = (x_.') * S * x_;

% define norm part
nor22 = epsilon * ( x_(1)^2 + x_(2)^2 );

% ----------- replace x_ with px ------------
% dJ = subs(dJ, x_(i), px(i));
% J = subs(J, x_(i), px(i));
% nor22 = subs(nor22, x_(i), px(i));

% define a sos program
Program1 = sosprogram(x_);

% define h(x_) as sums of squares
[Program1, h] = sossosvar(Program1, x_);

% add inequality constraint
Program1 = sosineq(Program1, -dJ - h*(rho - J) - nor22);

% set solver option

Program1 = sossolve(Program1);

SOLV = sosgetsol(Program1, h)

% OUTPUT IS REASONABLE
