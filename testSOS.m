%% example 1

clear; echo on;
syms x1 x2;

vartable = [x1, x2];

prog = sosprogram(vartable);

p = 2*x1^4 + 2*x1^3*x2 - x1^2*x2^2 + 5*x2^4;

prog = sosineq(prog, p);

prog = sossolve(prog);

echo off;

%% example 2
clear; echo on;
syms x1 x2 x3;
vars = [x1, x2, x3];

f = [-x1^3-x1*x3^2;
     -x2-x1^2*x2;
     -x3+3*x1^2*x3-3*x3/(x3^2+1)];
 
prog = sosprogram(vars);

[prog,V] = sospolyvar(prog,[x1^2; x2^2; x3^2],'wscoeff');

prog = sosineq(prog,V-(x1^2+x2^2+x3^2));

expr = -(diff(V,x1)*f(1)+diff(V,x2)*f(2)+diff(V,x3)*f(3))*(x3^2+1);
prog = sosineq(prog,expr);

prog = sossolve(prog);

SOLV = sosgetsol(prog,V)

echo off;