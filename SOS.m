clear; echo on;
pvar px1 px2 px3;
px = [px1; px2; px3];
mypi = 3.14159265358;

syms test;

x = sym('x', [3 1]);
u = sym('u', [2 1]);
x_ = sym('x_', [3 1]);
u_ = sym('u', [3 1]);
eq_x = [12.5; 12.5; mypi/4.0 ];
eq_u = [ 0.0; 0.0 ];

f_original = [ u(1) * cos(x(3)) ; u(1) * sin(x(3)) ; u(2)];

% tolerance 
epsilon = 0.001;
% initial value of rho
rho = 0.0001;

Q = [  10    0    0   ;
       0    10    0   ;
       0    0    5   ];
R = [  3    0   ;
       0    5   ];

[A_l, B_l] = sym_linearization(f_original, x, u, eq_x, eq_u);

[k, S, ~] = lqr(A_l, B_l, Q, R);

k
S

cg = cos(mypi/4.0);
sg = sin(mypi/4.0);

% define f_hat

% closed loop function
f_cl = toCL(f_original, x, u, eq_x, eq_u, k);

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
nor22 = epsilon * ( (x_(1))^2 + (x_(2))^2 + (x_(3))^2 );

% ----------- replace x_ with px ------------
% dJ = subs(dJ, x_(i), px(i));
% J = subs(J, x_(i), px(i));
% nor22 = subs(nor22, x_(i), px(i));

% define a sos program
Program1 = sosprogram(x_);

% define h(x_) as sums of squares
[Program1, h] = sossosvar(Program1, x_);

% add inequality constraint
Program1 = sosineq(Program1, (-dJ - h*(rho - J) - nor22) );

% call solver
Program1 = sossolve(Program1);





