function f = objective(z,u)

mode = 2;

if mode == 1
    R = 1.0;

    % Want slopes near zero, so minimize slope squared:
    f = 1.0 + u.^2*R;
    
elseif mode == 2
    R = 10.0;
    
    f = 1.0 + u.^2*R;
    
else
    f = 1.0 + u.^2;
end
    
end