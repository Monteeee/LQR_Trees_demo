function [a_z] = atomic_z(p, x)

order = zeros(length(x), 1);
terms = [];
for i = 1:length(x)
    order(i) = feval(symengine, 'degree', p, x(i));
    if order(i) > 0
        for j = 1:order(i)
           terms = [terms;x(i)];
        end
    end
end

num = length(terms);
if num == 1
    a_z = [terms(1), 1];
else
   % a_z = [prod(terms(1: fix(num/2))), prod(terms(fix(num/2) + 1:end)), 1];
   a_z = [terms(1), prod(terms(2:end)), 1];
end

end

