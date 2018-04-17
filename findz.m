load('poly_dJ.mat');
syms x_1 x_2 x_3 x_4;
x_bar = [x_1, x_2, x_3, x_4];
% dJ = 1 + x1*x2 + x1^2*x3 + x2*x3*x4^2 + x1^4 + x3^2*x4^2 + x3;
% dJ = x2*x3*x4^2;
tic
z = P_recursive(dJ, x_bar);
[c, t] = coeffs(dJ);
z
sos = z.' * z;
sos = sos(:);
all(ismember(t, sos))
toc
% for i = 1:length(z)
%     sos = z(2:end).' * z(2:end);
%     sos = sos(:);
%     if ( all(ismember(i, sos)) )
%         z = z(2:end);
%     else
%         z = [z(2:end), z(1)];
%     end
% end
% z
% sos = z.' * z;
% sos = sos(:);
% all(ismember(t, sos))