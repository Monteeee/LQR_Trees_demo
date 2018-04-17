function [monomials] = MinCoverSOS(poly1, max_order, x_)
% given a polynomial and a maximal acceptable order, this function tries to
% find a Sums of Squares to whichs covers all terms of poly1 while
% possessing least order (Note that this does not equals to least number of 
% terms, to get that a better method has to be made).
    
%     max_order = zeros(1, length(x_));
%     for i = 1:length(x_)
%         max_order(i) = feval(symengine, 'degree', poly, x_(i));
%     end

    l = length(max_order);
    if l == 1
        u1 = max_order;
        u2 = max_order;
        u3 = max_order;
        u4 = max_order;
    else
        u1 = max_order(1);
        u2 = max_order(2);
        u3 = max_order(3);
        u4 = max_order(4);
    end
        
    [~, terms_poly] = coeffs(poly1, x_);
    
    cell_x = num2cell(x_);
    monomials = [];
    flag = 0;
    for j_1 = 1:u1
        for j_2 = 1:u2
            for j_3 = 1:u3
                for j_4 = 1:u4
                    monomials = mpmonomials(cell_x, {0:j_1, 0:j_2, 0:j_3, 0:j_4});
                    sos = monomials * monomials.';
                    sos = sos(:);

                    if ( all(ismember(terms_poly, sos)) )
                        disp("found a feasible SOS, returning the monomials");
                        flag = 1;
                        break;
                    end
                end
                if (flag == 1)
                    break;
                end
            end
            if (flag == 1)
                    break;
            end
        end
        if (flag == 1)
                    break;
        end
    end
    
    for t = 1:length(monomials)
        sos = monomials(2:end) * monomials(2:end).';
        sos = sos(:);
        if ( all(ismember(terms_poly, sos)) )
            monomials = monomials(2:end);
        else
            monomials = [monomials(2:end);monomials(1)];
        end
    end
            
    % with the sos found, do some pruning, we can get the least terms
    % representation! But still it is very slow!
    
    disp("returning monomials");
    return;
end

