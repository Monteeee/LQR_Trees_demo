function [fcl_t] = symbolic_Taylor(fcl, x, xg, order)
%% *expired, use built-in function taylor instead. This function provides 2nd order taylor expansion


Nx = length(x);

x_bar = sym('x_', [Nx, 1]);

fcl_t = fcl;

for k=1:length(fcl)
    
    for i=1:Nx
        if i == 1
            firstpart = x_bar(1) * diff( fcl(k), x(1) );
        else
            firstpart = firstpart + x_bar(i) * diff( fcl(k), x(i) );
        end
    end
    
    for i=1:Nx
        for j=1:Nx
            if( i == 1 && j == 1)
                secondpart = 0.5 * x_bar(1) * x_bar(1) * ( diff(fcl(k), x(1), x(1)) );               
            else
                secondpart = secondpart + 0.5 * x_bar(i) * x_bar(j) * ( diff(fcl(k), x(i), x(j)) ); 
            end
        end
    end
    
    for i=1:Nx
        firstpart = subs(firstpart, x(i), xg(i));
        secondpart = subs(secondpart, x(i), xg(i));
    end
    
    fcl_t(k) = firstpart + secondpart;
end

fcl_t = vpa(fcl_t, 5);

end