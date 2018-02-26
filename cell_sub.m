function [ C3 ] = cell_sub( C1 , C2 )

n1 = size(C1, 1);
n2 = size(C1, 2);

if n1 ~= size(C2, 1) && n2 ~= size(C2, 2)
    error("sizes do not match!");
end
    
if isa(C1, 'cell')
    C3 = cell(size(C1));
    for i = 1:n1
        for j = 1:n2
            C3{i, j} = C1(i, j) - C2{i, j};
        end
    end
else
    C3 = zeros(size(C1));
    for i = 1:n1
        for j = 1:n2
            C3(i, j) = C1(i, j) - C2(i, j);
        end
    end
end
end