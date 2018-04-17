function [dSdt] = TV_Riccarti(t, S)

load('opti_traj.mat', 'soln', 'A_t', 'B_t', 'x_t', 'u_t', 'Q', 'R');

%% testing data loading and checking the trajectory

% t = linspace(soln.grid.time(1), soln.grid.time(end), 150);
% z_data = soln.interp.state(t);
% x_data = z_data(1,:);
% y_data = z_data(3,:);
% u_data = soln.interp.control(t);
% 
% plot(x_data,'DisplayName','x_data');hold on;plot(y_data,'DisplayName','y_data');hold off;


%% making Time-Variant Riccarti equation

fprintf("time: %d\n", t);

z_data = soln.interp.state(t);
x_data = z_data(1, :);
dx_data = z_data(2, :);
y_data = z_data(3, :);
dy_data = z_data(4, :);
u_data = soln.interp.control(t);

S = reshape(S, size(A_t)); %Convert from "n^2"-by-1 to "n"-by-"n"

% this will do the substitution much faster when x_data is vector
% f_At([x_t;u_t]) = A_t;
% val_A = f_At(x_data, dx_data, y_data, dy_data, u_data);
% 
% f_Bt([x_t;u_t]) = B_t;
% val_B = f_Bt(x_data, dx_data, y_data, dy_data, u_data);

val_A = double(subs(A_t, [x_t;u_t], [x_data;dx_data;y_data;dy_data;u_data]));
val_B = double(subs(B_t, [x_t;u_t], [x_data;dx_data;y_data;dy_data;u_data]));

% B_trans = val_B.';
% A_trans = val_A.';

dSdt = - (Q - S * ( 0.1 .* val_B * val_B.') * S + S*val_A + val_A.' * S);
% if size(u_data) == 1
%    Q_t = Q;
% else
%    Q_t = elem2vec(Q, size(u_data));
% end

% temp1 = cell_multi(val_B, B_trans);
% 
% temp2 = cell_multi(S, val_A);
% temp3 = cell_multi(A_trans, S);
% 
% temp4 = cell_multi(S, temp1);
% temp5 = cell_multi(temp4, S);
% 
% temp6 = cell_add(cell_add(temp5, temp2), temp3);
% 
% temp7 = cell_sub(Q_t, temp6);

% if isa(temp7, 'double')
%     dSdt = -1 .* temp7;
% end

dSdt = dSdt(:); %Convert from "n"-by-"n" to "n^2"-by-1

end