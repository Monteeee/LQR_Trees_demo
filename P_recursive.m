function [z] = P_recursive(p, x)
% UNTITLED3 Summary of this function goes here

[c, t] = coeffs(p);

nc = length(c);

if nc == 1
    z = atomic_z(p, x);
else
    % divide and conquer
    p1 = c(1: fix(nc/2)) * (t(1: fix(nc/2)).');
    z1 = P_recursive(p1, x);
    s1 = z1.' * z1;
    s1 = s1(:);
    
    if all(ismember(t(fix(nc/2) + 1:end), s1))
        z = z1;
        return;
    else
        p2 = c(fix(nc/2) + 1: end) * (t(fix(nc/2) + 1:end).');
        z2 = P_recursive(p2, x);
        z = union(z1, z2);
        return;
    end
end

