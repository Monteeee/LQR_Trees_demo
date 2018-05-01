clc; clear all;
% result from LTI LQR at the goal point

X0 = 1e3 * [0.0171    0.0268   -0.1279   -0.0409;
    0.0268    0.0786   -0.3968   -0.1270;
   -0.1279   -0.3968    5.4747    1.6305;
   -0.0409   -0.1270    1.6305    0.4888];

[T X] = ode45(@(t,X)TV_Riccarti(t, X), [9.85 0], X0);

Ns = sqrt(size(X, 2));

%% plot

pick = [1 2 3];

data1 = X(:, pick(1)) ./ (mean(X(:, pick(1))));
fprintf("scale of data 1: %d \n", mean(X(:, pick(1))));
data2 = X(:, pick(2)) ./ (mean(X(:, pick(2))));
fprintf("scale of data 2: %d \n", mean(X(:, pick(2))));
data3 = X(:, pick(3)) ./ (mean(X(:, pick(3))));
fprintf("scale of data 3: %d \n", mean(X(:, pick(3))));

% plot(T, data1, 'ro'); hold on;
% plot(T, data2, 'b*');
% plot(T, data3, 'g.-');

%% LTV LQR
load('opti_traj.mat', 'soln', 'B_t', 'x_t', 'u_t', 'A_t', 'Q', 'R');

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

%% preparation work for LTV verification
% get dS/dx = -Q + S*B*inv(R)*trans(B)*S + S*A + trans(A)*S
dSdt = zeros(Ns, Ns, length(t_temp));
for i=1:length(t_temp)
    Stemp = reshape(X(i, :), [4 4]);
    dSdt(:,:,i) = -Q + Stemp * B_mat(i, :).' * 0.1 * B_mat(i, :) * Stemp + ...
                    Stemp * A_mat(:, :, i) + A_mat(:, :, i).' * Stemp;
end

d1 = sqrt(size(X, 2));
d2 = d1;
% use cubic spline to interpolate the S matrix elements

% Xt = cell(d1, 1);
% 
% for i = 1:d1
%     pp = spline(T, z_data(i, :));
%     Xt{i} = pp.coefs;
% end
% 
% pp = spline(T, u_data);
% Ut = pp.coefs;
St = cell(d1, d2);
dSt = cell(d1, d2);
for i = 1:d1
    for j = 1:d2
        pp = spline(T, X(:, j + (i-1)*4));
        St{i, j} = pp.coefs;
        
        pp = spline(T, dSdt(i, j, :));
        dSt{i, j} = pp.coefs;
    end
end

% spline coefficients:
% a(t−t1)^3 + b(t−t1)^2 + c(t−t1) + d, t1 is the starting point of the
% time interval

%% use spline-like S elements to get J_hat
x_ = sym('x_', [d1, 1]);
syms t;
J_hat_t = sym(zeros(length(T)-1, 1));
for tk = 1:length(T)-1
    for i = 1:d1
        for j = 1:d2
            coff_t = St{i, j}(tk, :);
            t1 = T(tk);
            sst = coff_t(1)*(t - t1)^3 + coff_t(2)*(t - t1)^2 + coff_t(3)*(t - t1) + coff_t(4);
            J_hat_t(tk) = J_hat_t(tk) + ( x_(i) * x_(j) * sst );
        end
    end
end


%% Get Taylor expansion of f on each Time knot
x = sym('x', [4 1]);
u = sym('u', [1 1]);

m1 = 2.0;
m2 = 0.5;
g = 9.81;
l = 1.0;
d = 1.0;

f_original = [ x(2) ;
               ( l*m2*sin(x(3))*x(4)^2 + u(1) + m2*g*cos(x(3))*sin(x(3)) ) / ( m1 + m2*(1 - (cos(x(3)))^2) );
               x(4) ;
               - ( l*m2*cos(x(3))*sin(x(3))*x(4)^2 + u(1)*cos(x(3)) + (m1+m2)*g*sin(x(3)) ) / ( l*m1 + l*m2*(1 - (cos(x(3)))^2) ) ];
                    
f_t = sym(zeros(length(f_original), length(T)));

% taylor expansion order
N_taylor = 2;

for tk = 1:length(T)
    eq_x = [x_data(tk);dx_data(tk);y_data(tk);dy_data(tk)];
    output = taylor(f_original, x.', eq_x.', 'order', N_taylor);
    % substitute to get the closed loop f
%     x_sub = sym(zeros(4, 1));
%     for i = 1:length(x)
%         coff_x = Xt{i}(tk, :);
%         t1 = T(tk);
%         x_sub(i) = x_(i) + coff_t(1)*(t - t1)^3 + coff_t(2)*(t - t1)^2 + coff_t(3)*(t - t1) + coff_t(4);
%     end
%     coff_u = Ut(tk, :);
%     u_sub = x_(i) + coff_t(1)*(t - t1)^3 + coff_t(2)*(t - t1)^2 + coff_t(3)*(t - t1) + coff_t(4)
    x_sub = z_data(:, tk) + x_;
    u_sub = u_data(tk) - K_t(tk, :) * x_;
    f_t(:, tk) = subs(output, [x;u], [x_sub;u_sub]);
end

%% get coefficients series from Taylor expansion of f (slow, because of symbolic)

% ###### what to do with so many coefficients to get spline? #####

f_coeff_all = cell(4, 1);
terms_all = cell(4, 1);

for i = 1:size(f_t, 1)
    last_term = [];
    
    for tk = 1:length(T)
        fprintf("iteration %d \n", tk);
        
        if tk == 1
            [coeff_temp, terms] = coeffs(f_t(i, tk), x_);
            coefficients = zeros(length(terms), length(T));
            last_term = terms;
        else
            [coeff_temp, terms] = coeffs(f_t(i, tk), x_);
        end
        % disp(length(terms));
        
        if ~isequal(terms, last_term)
            disp("not equal")
            if tk == length(T)
                if isequal(terms, last_term(1:length(last_term)-1))
                    terms = [terms, sym(1)];
                    coeff_temp = [coeff_temp, sym(0)];
                end
                % terms
                % last_term
                % coeff_temp
            end
        end
        
        coefficients(:, tk) = coeff_temp(:);

%         if ~isequal(terms, last_term)
%             disp("different content!");
%             if tk ~= 1
%                 disp(terms == last_term);
% 
%             end
%         end

    end
    
    % keep copy of the terms we have for the coefficients
    terms_all{i} = last_term;
    f_coeff_all{i} = coefficients;
end
    
%% get spline of each coefficient over each time interval (doing spline is actually very fast)
coeff_ft = cell(size(f_coeff_all, 1), 1);
for j = 1:size(f_coeff_all, 1)
    coeff_fi = zeros(size(f_coeff_all{j}, 1), length(T)-1, 4);
    for i = 1:size(f_coeff_all{j}, 1)
        pp = spline(T, f_coeff_all{j}(i, :));
        coeff_fi(i, :, :) = pp.coefs;
    end
    coeff_ft{j} = coeff_fi;
end

ft_spline = cell(length(T)-1, size(f_coeff_all, 1));
for k = 1:size(f_coeff_all, 1) % 4
    for i = 1:length(T) - 1 % 96
        item = sym(0);
        for j = 1:length(terms_all{k})  % change with f function
            cft = coeff_ft{k}(j, i, :);
            t1 = T(i);
            ct = cft(1)*(t - t1)^3 + cft(2)*(t - t1)^2 + cft(3)*(t - t1) + cft(4);
            item = item + ct * terms_all{k}(j);
        end
        ft_spline{i, k} = item;
    end
end

%% get dJ_hat using all kinds of splines (abandoned due to ridiculous complexity)
dJ_hat_t2 = sym(zeros(length(T)-1, 1));
for tk = 1:length(T)-1
    for i = 1:d1
        for j = 1:d2
            coeff_St = St{i, j}(tk, :);
            coeff_dSt = dSt{i, j}(tk, :);
            t1 = T(tk);
            sst = coeff_St(1)*(t - t1)^3 + coeff_St(2)*(t - t1)^2 + coeff_St(3)*(t - t1) + coeff_St(4);
            dsst = coeff_dSt(1)*(t - t1)^3 + coeff_dSt(2)*(t - t1)^2 + coeff_dSt(3)*(t - t1) + coeff_dSt(4);
            
            dJ_hat_t2(tk) = dJ_hat_t2(tk) + ( x_(i) * x_(j) * dsst ) + 2 * x_(i) * sst *  ft_spline{tk, j};
        end
    end
end

%% get dJ_hat for each time instant first and then compute spline for it

dJ_hat_t = sym(zeros(length(T), 1));

for tk = 1:length(T)
   for i = 1:d1
       for j = 1:d2
            dJ_hat_t(tk) = dJ_hat_t(tk) + x_(i) * x_(j) * dSdt(i, j, tk) + ...
                2 * x_(i) * X(tk, j + (i-1)*4) * f_t(j, tk);
       end
   end
   if tk == 1
        [~, terms_of_dJ] = coeffs(dJ_hat_t(tk), x_);
   end
   % following is a check for terms of dJ_hat_t, if all of them are the
   % same then we can safely do the spline
%    [coeff_temp, terms] = coeffs(dJ_hat_t(tk), x_);
%    if tk == 1
%        last_terms = terms;
%    elseif ~isequal(last_terms, terms)
%        disp("not equal terms of dJ!");
%    end
end

coeff_dJ = zeros(length(T), length(terms_of_dJ));
for tk = 1:length(T)
   [coeff_dJ(tk, :), ~] = coeffs(dJ_hat_t(tk), x_);
end

sp_dJ = zeros(length(terms_of_dJ), length(T)-1, 4);
for i = 1:length(terms_of_dJ)
    pp = spline( T, coeff_dJ(:, i) );
    sp_dJ(i, :, :) = pp.coefs;
end

dJ_spline = cell(length(T)-1, 1);
for i = 1:length(T) - 1 % 96
    item = sym(0);
    for j = 1:length(terms_of_dJ)
        cft = sp_dJ(j, i, :);
        t1 = T(i);
        ct = cft(1)*(t - t1)^3 + cft(2)*(t - t1)^2 + cft(3)*(t - t1) + cft(4);
        item = item + ct * terms_of_dJ(j);
    end
    dJ_spline{i} = item;
end

%% the SOS programming for (39) in Tedrake paper 2010
beta1k = 0.1;
beta0k = 0.1;
rho_k = beta1k * t + beta0k;
drho_k = beta1k;

dJ_hat_trunc = dJ_spline;

k = 30;
T_r = fliplr(T);

[c1, term1] = coeffs(dJ_spline{k}, [x_;t]);
c1 = double(c1);

% maybe try to do truncatoin here !
c1_trunc = c1(abs(c1) > 0.01 * mean(abs(c1)));
term1_trunc = term1(abs(c1) > 0.01 * mean(abs(c1)));

dJ_hat_trunc{k} = c1_trunc * term1_trunc.';

%% 
% define a sos program
Program1 = sosprogram([x_;t]);
tic
mz = P_recursive(dJ_hat_trunc{k}, [x_;t]);

% define h1(x_, t), h2(x_, t), h3(x_, t) as sums of squares

[Program1, h2] = sossosvar(Program1, mz.');
[Program1, h3] = sossosvar(Program1, mz.');

% VEC_x = monomials(x_, [1 2 3 4]);
% VEC_t = monomials(t, [1 2 3 4 5 6]);
% VEC_ = VEC_x * VEC_t.';
% VEC_ = reshape(VEC_, [numel(VEC_), 1]);

% disp("adding h1");
[Program1, h1] = sospolyvar(Program1, term1_trunc.');
toc
%% constraint
disp("adding constraint");
tic
% add inequality constraint
Program1 = sosineq(Program1, -(dJ_hat_trunc{k} - beta1k + h1 * (rho_k - J_hat_t(k)) + ...
                             h2 * (t - T_r(k)) + h3 * (T_r(k+1) - t) ));
toc
tic
% set solver option
disp("start solving");
Program1 = sossolve(Program1);

SOLV1 = sosgetsol(Program1, h1);
disp(SOLV1)
SOLV2 = sosgetsol(Program1, h2);
disp(SOLV2)
SOLV3 = sosgetsol(Program1, h3);
disp(SOLV3)
toc

%% 
disp("done!");

