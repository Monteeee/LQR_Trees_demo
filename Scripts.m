A = zeros(3);
sq2 = sqrt(2);
B = [sq2 0; sq2 0; 0 1];
Q = [  10    0    0   ;
       0    10    0   ;
       0    0    5   ];
R = [  3    0   ;
       0    5   ];

[K,S,e] = lqr(A,B,Q,R);

% disp(K);
% disp(S);
% disp(e);

mypi = 3.14159265358;

x = sym('x', [3 1]);
u = sym('u', [2 1]);

f_original = [ u(1) * cos(x(3)) ; u(1) * sin(x(3)) ; u(2)];

eq_x = [12.5; 12.5; mypi/4.0 ];
eq_u = [ 0.0; 0.0 ];

[A_l, B_l] = sym_linearization(f_original, x, u, eq_x, eq_u);

%B_l

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