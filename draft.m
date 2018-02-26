% result from LTI LQR at the goal point

X0 = 1e3 * [0.0171    0.0268   -0.1279   -0.0409;
    0.0268    0.0786   -0.3968   -0.1270;
   -0.1279   -0.3968    5.4747    1.6305;
   -0.0409   -0.1270    1.6305    0.4888];

[T X] = ode45(@(t,X)TV_Riccarti(t, X), [9.85 0], X0);

%% plot

pick = [1 2 3];

data1 = X(:, pick(1)) ./ (mean(X(:, pick(1))));
fprintf("scale of data 1: %d \n", mean(X(:, pick(1))));
data2 = X(:, pick(2)) ./ (mean(X(:, pick(2))));
fprintf("scale of data 2: %d \n", mean(X(:, pick(2))));
data3 = X(:, pick(3)) ./ (mean(X(:, pick(3))));
fprintf("scale of data 3: %d \n", mean(X(:, pick(3))));

plot(T, data1, 'ro'); hold on;
plot(T, data2, 'b*');
plot(T, data3, 'g.-');

%% LTV LQR
load('opti_traj.mat', 'soln', 'B_t', 'x_t', 'u_t');

% get B(t), we need it to compute K(t)
t = T.';
z_data = soln.interp.state(t);
x_data = z_data(1, :);
dx_data = z_data(2, :);
y_data = z_data(3, :);
dy_data = z_data(4, :);
u_data = soln.interp.control(t);

f_Bt([x_t;u_t]) = B_t;
val_B = f_Bt(x_data, dx_data, y_data, dy_data, u_data);

K_t = zeros(length(t), 4);

B_mat = zeros(length(t), 4);
for i=1:length(t)
    B_mat(i, 1) = double(val_B{1}(i));
    B_mat(i, 2) = double(val_B{2}(i));
    B_mat(i, 3) = double(val_B{3}(i));
    B_mat(i, 4) = double(val_B{4}(i));
end

for i=1:length(t)
    K_t(i, :) = 0.1 .* (B_mat(i, :) * reshape(X(i, :), [4 4]));
end

disp("done!");