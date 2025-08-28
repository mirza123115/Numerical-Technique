clear all
close all
clc

% Exercise 1: Polynomial operations
% Define A(x) = 3x^2 + 2x - 4
A = [3 2 -4];

% Define B(x) = 2x^3 - 2
B = [2 0 0 -2];

% Multiply polynomials C(x) = A(x)*B(x)
C = conv(A, B);

disp('Polynomial A(x):');
disp(A);
disp('Polynomial B(x):');
disp(B);
disp('Polynomial C(x) = A(x)*B(x):');
disp(C);

% Find roots of A(x), B(x), C(x)
roots_A = roots(A);
roots_B = roots(B);
roots_C = roots(C);

disp('Roots of A(x):');
disp(roots_A);
disp('Roots of B(x):');
disp(roots_B);
disp('Roots of C(x):');
disp(roots_C);


%% Data from Table 1
x = [0 1 2 3 4 5 6];
f = [0 0.8415 0.9093 0.1411 -0.7568 -0.9589 -0.2794];

% Create figure
figure; hold on; grid on;

% Plot original data points
plot(x, f, 'ro', 'MarkerFaceColor', 'r');

% Piecewise linear interpolation using 2-point formula
for i = 1:length(x)-1
    % Take two consecutive points
    x0 = x(i);   f0 = f(i);
    x1 = x(i+1); f1 = f(i+1);
    
    % Generate fine points between x0 and x1
    xx = linspace(x0, x1, 50);  
    
    % Linear interpolation formula
    yy = f0 + ( (f1 - f0) / (x1 - x0) ) * (xx - x0);
    
    % Plot line segment
    plot(xx, yy, 'b-');
end

title('Linear Interpolation (Piecewise)');
xlabel('x');
ylabel('f(x)');
legend('Data Points', 'Interpolated Lines');


%% Exercise 3: Interpolation of sin(x)

% Original data points
x = 0:10;           
y = sin(x);         

% Points where we want interpolated values
xi = 0:0.25:10;     

% Linear interpolation
yi = interp1(x, y, xi, 'linear');

% Plot results
figure; hold on; grid on;

% Original points
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Data Points');

% Interpolated curve
plot(xi, yi, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Linear Interpolant');

% Actual sine curve for comparison
xfine = 0:0.01:10;
yfine = sin(xfine);
plot(xfine, yfine, 'g--', 'LineWidth', 1, 'DisplayName', 'Actual sin(x)');

title('Exercise 3: Linear Interpolation of sin(x)');
xlabel('x'); ylabel('y');
legend show;







%% Exercise 4 & 5: Lagrange Polynomial Interpolation
clear; clc; close all;

%% Step 1: Define data points (example data)
x = [0 1 2 3 4 5 6];          % x-coordinates
f = [0 0.9038 0.2255 -0.3577 0.07321 -0.00313 -0.0001521];  % f(x) values

N = length(x) - 1;            % Degree of polynomial

%% Step 2: Define interpolation point
x_interp = 2.5;               % You can change this value

%% Step 3: Lagrange Polynomial Algorithm
SUM = 0;                       % Initialize sum

for i = 0:N
    P = 1;                     % Initialize product
    for j = 0:N
        if j ~= i
            P = P * (x_interp - x(j+1)) / (x(i+1) - x(j+1)); % Lagrange basis
        end
    end
    SUM = SUM + f(i+1) * P;     % Add weighted basis to sum
end

fprintf('Interpolated value at x = %.2f is %.4f\n', x_interp, SUM);

%% Step 4: Plot Lagrange interpolating polynomial
x_plot = linspace(min(x), max(x), 200);
y_plot = zeros(size(x_plot));

for k = 1:length(x_plot)
    xi = x_plot(k);
    SUM_k = 0;
    for i = 0:N
        P = 1;
        for j = 0:N
            if j ~= i
                P = P * (xi - x(j+1)) / (x(i+1) - x(j+1));
            end
        end
        SUM_k = SUM_k + f(i+1) * P;
    end
    y_plot(k) = SUM_k;
end

figure;
plot(x_plot, y_plot, 'b-', 'LineWidth', 2); hold on;
plot(x, f, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('x'); ylabel('f(x)');
title('Lagrange Interpolating Polynomial');
grid on;
legend('Interpolating Polynomial','Data Points');





%% Exercise 6 & 7: Newton Polynomial Interpolation
clear; clc; close all;

%% Step 1: Define data points (example from table 1)
x = [0 1 2 3 4 5 6];             % x-coordinates
f = [0 0.9038 0.2255 -0.3577 0.07321 -0.00313 -0.0001521];  % f(x) values

N = length(x);                    % Number of points

%% Step 2: Compute Divided Differences Table
D = zeros(N,N);
D(:,1) = f(:);                    % First column is f(x)

for j = 2:N
    for i = j:N
        D(i,j) = (D(i,j-1) - D(i-1,j-1)) / (x(i) - x(i-j+1));
    end
end

%% Step 3: Evaluate Newton Polynomial at a given point
x_interp = 2.5;                   % Interpolation point
P = D(1,1);                        % Initialize with first term
prod_term = 1;

for j = 2:N
    prod_term = prod_term * (x_interp - x(j-1));
    P = P + D(j,j) * prod_term;
end

fprintf('Interpolated value at x = %.2f is %.4f\n', x_interp, P);

%% Step 4: Plot Newton interpolating polynomial
x_plot = linspace(min(x), max(x), 200);
y_plot = zeros(size(x_plot));

for k = 1:length(x_plot)
    xi = x_plot(k);
    Pk = D(1,1);
    prod_term = 1;
    for j = 2:N
        prod_term = prod_term * (xi - x(j-1));
        Pk = Pk + D(j,j) * prod_term;
    end
    y_plot(k) = Pk;
end

figure;
plot(x_plot, y_plot, 'b-', 'LineWidth', 2); hold on;
plot(x, f, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('x'); ylabel('f(x)');
title('Newton Interpolating Polynomial');
grid on;
legend('Newton Polynomial','Data Points');



