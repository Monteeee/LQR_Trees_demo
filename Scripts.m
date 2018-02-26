%% basic setting of problem
A = zeros(3);
sq2 = sqrt(2);
B = [sq2 0; sq2 0; 0 1];
Q = [  10    0    0   0   0;
       0    10    0   0   0;
       0     0    10   0   0;
       0     0    0   10   0;
       0     0    0   0   10];
R = [  10    0   ;
       0    10   ];

% disp(K);
% disp(S);
% disp(e);

mypi = 3.14159265358;

x = sym('x', [5 1]);
u = sym('u', [2 1]);

r = 0.1;
l = 0.3;
I = 1.0;

f_original = [ r/2.0 * (x(4) + x(5)) * cos(x(3)) ; 
               r/2.0 * (x(4) + x(5)) * sin(x(3)) ;
               r/l * (x(4) - x(5));
               u(1)/I;
               u(2)/I];

eq_x = [12.5; 12.5; mypi/4.0; 0; 0 ];
eq_u = [ 0.0; 0.0 ];

%% linearization and more
[A_l, B_l] = sym_linearization(f_original, x, u, eq_x, eq_u);

A_l
B_l

[K_new, S_new, ~] = lqr(A_l, B_l, Q, R);

K_new
S_new

f_cl = toCL(f_original, x, u, eq_x, eq_u, K_new);

f_cl_t = symbolic_Taylor(f_cl, x, eq_x);

for i=1:length(f_cl_t)
    [C, T] = coeffs(f_cl_t(i));
    C(abs(C) < 1e-7) = 0;
    f_cl_t(i) = dot(C, T);
end
f_cl_t = vpa(f_cl_t, 5);
f_cl_t