function [M3] = cell_multi(M1,M2)

n11 = size(M1, 1);
n12 = size(M1, 2);

n21 = size(M2, 1);
n22 = size(M2, 2);

if n12 ~= n21
    error("the dimenstions do not match!");
end

if isa(M1, 'cell')
    M3 = cell(n11, n22);

    sum = zeros(size(M1{1, 1}));

    for i=1:n11
        for j=1:n22
            sum = zeros(size(M1{1, 1}));
            for k=1:n12
                sum = sum + M1{i, k} .* M2{k, j};
            end
            M3{i, j} = sum;
        end
    end
    
else
    M3 = zeros(n11, n22);

    for i=1:n11
        for j=1:n22
            sum = 0;
            for k=1:n12
                sum = sum + M1(i, k) .* M2(k, j);
            end
            M3(i, j) = sum;
        end
    end
end

end

