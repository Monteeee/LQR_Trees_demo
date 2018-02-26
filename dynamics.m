function dz = dynamics(z,u)

mode = 2;

if mode == 1
    
    x = z(1, :);
    v = z(2, :);
    k = 0.5;
    m = 1.0;

    dz = [v; -k/m * x + 1/m * u];

elseif mode == 2
    
    m1 = 2.0;
    m2 = 0.5;
    g = 9.81;
    l = 1.0;
    d = 1.0;
    
    x = z;
    
    dz = [ x(2, :) ;
         ( l*m2*sin(x(3, :)).*x(4, :).^2 + u + m2*g*cos(x(3, :)).*sin(x(3, :)) ) / ( m1 + m2*(1 - (cos(x(3))).^2) );
           x(4, :) ;
          -( l*m2*cos(x(3, :)).*sin(x(3, :)).*x(4).^2 + u.*cos(x(3, :)) + (m1+m2)*g.*sin(x(3, :)) ) / ( l*m1 + l*m2.*(1 - (cos(x(3))).^2) ) ];  
else
    dz = [];
end
    
    
end