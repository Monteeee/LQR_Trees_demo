function [r, alpha, beta] = xytopolar(dx, dy, theta)

r = sqrt(dx^2 + dy^2);

alpha = -theta + atan(dy, dx);

alpha = mod(alpha, 3.14159265358 / 2.0);

beta = -theta - alpha;

beta = mod(beta, 3.14159265358);

end

