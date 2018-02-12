function [fcl] = toCL(f, x, u, xg, ug, K)
%% provide closed-loop function

Nu = length(ug);

newu = ug - K*(x - xg);

fcl = f;

for i=1:Nu
   fcl = subs(fcl, u(i), newu(i));
end

fcl = vpa(fcl, 5);

end
