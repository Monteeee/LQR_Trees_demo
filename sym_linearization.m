function [A, B] = sym_linearization(f, x, u, eq_x, eq_u)

Nx = length(f);
Nu = length(eq_u);

A = sym('trash', [Nx Nx]);
B = sym('trash', [Nx Nu]);

for i=1:Nx
    for j=1:Nx
        A(i, j) = diff(f(i), x(j));
    end
end

for i=1:Nx
    for j=1:Nu
        B(i, j) = diff(f(i), u(j));
    end
end

A = subs(A, x, eq_x);
A = subs(A, u, eq_u);
A = double(A);
B = subs(B, x, eq_x);
B = subs(B, u, eq_u);
B = double(B);
end

