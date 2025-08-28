%% Linear Least Squares Regression - Table 1
clc; clear; close all;

%% Step 1: Given data (from Table 1)
x = [1 2 3 4 5 6 7];
y = [0.5 2.5 2.0 4.0 3.5 6.0 5.5];

n = length(x);  % number of data points

%% Step 2: Compute sums for coefficients
sumx = sum(x);
sumy = sum(y);
sumxy = sum(x .* y);
sumxsq = sum(x .^2);

%% Step 3: Compute linear regression coefficients
a1 = (n*sumxy - sumx*sumy) / (n*sumxsq - sumx^2);  % slope
a0 = mean(y) - a1 * mean(x);                        % intercept

fprintf('Linear regression coefficients:\n');
fprintf('a0 (intercept) = %.8f\n', a0);
fprintf('a1 (slope) = %.8f\n', a1);

%% Step 4: Plot data points and fitted line
figure;
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on;
xlabel('x'); ylabel('y');
title('Linear Least Squares Regression - Table 1');
grid on;

% Fitted line
ym = a0 + a1 * x;
plot(x, ym, 'b-', 'LineWidth', 2);
legend('Observed Data','Fitted Line');






%% Polynomial Regression (Second Order) - Table 2
clc; clear; close all;

%% Step 1: Given data (Table 2)
x = [0 1 2 3 4 5];
y = [2.1 7.7 13.6 27.2 40.9 61.1];

n = length(x);

%% Step 2: Construct matrix for normal equations
% Sum calculations
sumx0 = n;              % sum(x^0)
sumx1 = sum(x);          % sum(x)
sumx2 = sum(x.^2);       % sum(x^2)
sumx3 = sum(x.^3);       % sum(x^3)
sumx4 = sum(x.^4);       % sum(x^4)

sumy0 = sum(y);           % sum(y)
sumxy = sum(x .* y);      % sum(x*y)
sumx2y = sum((x.^2) .* y);% sum(x^2*y)

% Normal equations in matrix form
A = [sumx0 sumx1 sumx2;
     sumx1 sumx2 sumx3;
     sumx2 sumx3 sumx4];

B = [sumy0; sumxy; sumx2y];

%% Step 3: Solve for coefficients using matrix inversion
coeff = A\B;  % coeff = [a0; a1; a2]

a0 = coeff(1);
a1 = coeff(2);
a2 = coeff(3);

fprintf('Polynomial coefficients:\n');
fprintf('a0 = %.5f\n', a0);
fprintf('a1 = %.5f\n', a1);
fprintf('a2 = %.5f\n', a2);

%% Step 4: Plot data points and fitted polynomial
x_plot = linspace(min(x), max(x), 200);
y_plot = a0 + a1*x_plot + a2*x_plot.^2;

figure;
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on;
plot(x_plot, y_plot, 'b-', 'LineWidth', 2);
xlabel('x'); ylabel('y');
title('Second Order Polynomial Regression - Table 2');
grid on;
legend('Data Points','Fitted Polynomial');






%% Linearization of Nonlinear Power Law Model
clc; clear; close all;

%% Step 1: Given data (Table 3)
x = [1 2 3 4 5];
y = [0.5 1.7 3.4 5.7 8.4];

%% Step 2: Transform data using base 10 logarithm
logx = log10(x);
logy = log10(y);

%% Step 3: Perform linear regression on log-log data
n = length(x);

sum_logx = sum(logx);
sum_logy = sum(logy);
sum_logx_logy = sum(logx .* logy);
sum_logx2 = sum(logx.^2);

% Slope (b2) and intercept (loga2)
b2 = (n*sum_logx_logy - sum_logx*sum_logy) / (n*sum_logx2 - sum_logx^2);
loga2 = mean(logy) - b2 * mean(logx);

% Compute a2
a2 = 10^loga2;

fprintf('Power Law Model y = a2 * x^b2\n');
fprintf('a2 = %.4f\n', a2);
fprintf('b2 = %.4f\n', b2);

%% Step 4: Plot original data and fitted curve
x_fit = linspace(min(x), max(x), 100);
y_fit = a2 * x_fit.^b2;

figure;
plot(x, y, 'ro', 'MarkerFaceColor','r','MarkerSize',8); hold on;
plot(x_fit, y_fit, 'b-', 'LineWidth',2);
xlabel('x'); ylabel('y');
title('Power Law Fit using Linearization');
grid on;
legend('Data Points', 'Fitted Curve');
