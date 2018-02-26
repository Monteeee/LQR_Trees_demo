function [c_a, c_b] = finalconstraint(x0, xf)

mode = 2;

if mode == 1

    tree_point = [0.1, 0.0; 0.3, -0.05; -0.15, 0.08]';

    c_a = norm(xf-tree_point(:, 1)) * norm(xf-tree_point(:, 2)) * norm(xf-tree_point(:, 3));

    c_b = [];
    
elseif mode == 2
    d = 1.0;
    my_pi = 3.14159265358;
    tree_point = [d,        0.0, my_pi,        0.0;
                  d - 0.05, 0.0, my_pi - 0.05, 0.0;
                  d - 0.1, 0.08, my_pi + 0.03, -0.01 ]';
    
    mul = 1.0;
    for i = 1:size(tree_point, 2)
        mul = mul * norm(xf-tree_point(:, i));
    end
    
    c_a = mul;
    c_b = [];
else
    c_a = [];
    c_b = [];
end
              
end

