clc; clear all;
% result from LTI LQR at the goal point
load('opti_traj_1.mat', 'soln', 'B_t', 'x_t', 'u_t', 'A_t', 'Q', 'R', 'S');

X0 = S;

[T, X] = ode45(@(t,X)TV_Riccarti(t, X), [9.85 0], X0);

Ns = sqrt(size(X, 2));

%% LTV LQR

% get B(t), we need it to compute K(t)
t_temp = T.';
z_data = soln.interp.state(t_temp);
x_data = z_data(1, :);
dx_data = z_data(2, :);
y_data = z_data(3, :);
dy_data = z_data(4, :);
u_data = soln.interp.control(t_temp);

f_Bt([x_t;u_t]) = B_t;
val_B = f_Bt(x_data, dx_data, y_data, dy_data, u_data);

f_At([x_t;u_t]) = A_t;
val_A = f_At(x_data, dx_data, y_data, dy_data, u_data);

B_mat = zeros(length(t_temp), 4);
for i=1:length(t_temp)
    B_mat(i, 1) = double(val_B{1}(i));
    B_mat(i, 2) = double(val_B{2}(i));
    B_mat(i, 3) = double(val_B{3}(i));
    B_mat(i, 4) = double(val_B{4}(i));
end

A_mat = zeros(size(val_A, 1), size(val_A, 2), length(t_temp));
for i=1:length(t_temp)
    for j = 1:size(val_A, 1)
        for k = 1:size(val_A, 2)
            A_mat(j, k, i) = double(val_A{j, k}(i));
        end
    end
end

K_t = zeros(length(t_temp), 4);

for i=1:length(t_temp)
    K_t(i, :) = 0.1 .* (B_mat(i, :) * reshape(X(i, :), [4 4]));
end

%% run simulation

m1 = 2.0;
m2 = 0.5;
g = 9.81;
l = 1.0;
d = 1.0;
equilibrium_x = [d ; 0.0; pi; 0.0];
equilibrium_u = [ 0.0 ];
tree_point = [d,        0.0, pi,        0.0;
              d - 0.05, 0.0, pi - 0.05, 0.0;
              d - 0.1, 0.08, pi + 0.03, -0.01 ]';

T_inc = fliplr(T.');

K_t_data.time = fliplr(T.').';
K_t_data.signals.values = fliplr(K_t.').';
K_t_data.signals.dimensions = size(K_t, 2);

X_t_data.time = fliplr(T.').';
X_t_data.signals.values = fliplr(z_data).';
X_t_data.signals.dimensions = size(z_data, 1);

U_t_data.time = fliplr(T.').';
U_t_data.signals.values = fliplr(u_data).';
U_t_data.signals.dimensions = size(u_data, 1);

x0 = equilibrium_x(1) - 0.4;
dx0 = 0;
theta0 = equilibrium_x(3) - 0.4;
dtheta0 = 0;

simOut = sim('ddx');
states_out = states.Data;
time_out = states.Time;

t_end = T_inc(end);
[time_diff, idx] = min(abs(time_out - t_end));
states_diff = 100;
err_weight = [0.75; 0.75; 1.5; 1.5];
for i = 1:size(tree_point, 2)
    err_temp = norm( err_weight .* ...
        (states_out(1:4, :, idx) - tree_point(:, i)) );
    if states_diff > err_temp
        states_diff = err_temp;
        closest_point = tree_point(:, i);
    end
end

state_tolerance = 0.2;
if states_diff < state_tolerance
    fprintf("states final norm error: %f \n", states_diff);
    fprintf("error by state: \n");
    fprintf("- x     : %f \n", states_out(1, :, idx) - closest_point(1) );
    fprintf("- dx    : %f \n", states_out(2, :, idx) - closest_point(2) );
    fprintf("- theta : %f \n", states_out(3, :, idx) - closest_point(3) );
    fprintf("- dtheta: %f \n", states_out(4, :, idx) - closest_point(4) );
    disp('state is ok');
end

if max(abs(states_out(5, :, :))) < 50
    fprintf("max control signal: %f \n", max(abs(states_out(5, :, :))));
    disp('control signal is ok');
end

    