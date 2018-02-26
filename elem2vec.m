function [ext_M] = elem2vec(M, vec_size)

n1 = size(M, 1);
n2 = size(M, 2);


ext_M = cell(n1, n2);

for i = 1:n1
    for j = 1:n2
        ext_M{i, j} = M(i, j) .* ones(vec_size);
    end
end

end

