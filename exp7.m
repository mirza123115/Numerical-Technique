clc; clear;

% Tabulated data from Table 7.1
x = [1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8];
f = [4.953 6.050 7.389 9.025 11.023 13.468 16.445 20.086 24.533 29.964 36.598 44.701];

%% Composite Trapezoidal Rule
fprintf('Composite Trapezoidal Rule:\n');
h_values = [0.2 0.4 0.6];  % different step sizes

for h = h_values
    n = (x(end) - x(1)) / h;  % number of subintervals
    xi = x(1):h:x(end);
    
    % Interpolate f values for given h if necessary
    fi = interp1(x, f, xi, 'linear');  
    
    % Trapezoidal formula
    integral_trap = fi(1)/2 + fi(end)/2 + sum(fi(2:end-1));
    integral_trap = integral_trap * h;
    
    fprintf('h = %.1f -> Integral ≈ %.5f\n', h, integral_trap);
end

%% Simpson's 1/3 Rule (requires odd number of points)
fprintf('\nSimpson 1/3 Rule:\n');
h = 0.2;  % choose a step that gives odd number of points
n = (x(end-1) - x(1))/h;  % exclude last point to have odd points
xi = x(1):h:x(end-1);
fi = interp1(x, f, xi, 'linear');

if mod(length(fi),2)==1  % check for odd number of points
    integral_simp13 = fi(1) + fi(end) + 4*sum(fi(2:2:end-1)) + 2*sum(fi(3:2:end-2));
    integral_simp13 = integral_simp13 * h / 3;
    fprintf('Integral ≈ %.5f\n', integral_simp13);
else
    fprintf('Number of points not suitable for Simpson 1/3\n');
end


%% Simpson's 3/8 Rule (requires n divisible by 3)
fprintf('\nSimpson 3/8 Rule:\n');
h = 0.2;  % step for 3/8
n = (x(end-2) - x(1))/h;  % make n divisible by 3
xi = x(1):h:x(end-2);
fi = interp1(x, f, xi, 'linear');

if mod(length(fi)-1,3)==0
    integral_simp38 = fi(1) + fi(end) + 3*sum(fi([2:3:end-1 3:3:end-2])) + 2*sum(fi(4:3:end-3));
    integral_simp38 = integral_simp38 * 3*h / 8;
    fprintf('Integral ≈ %.5f\n', integral_simp38);
else
    fprintf('Number of points not suitable for Simpson 3/8\n');
end




%% Adaptive Trapezoidal Integration

clc; clear;

% Function to integrate
f = @(x) x.^2 .* exp(-x.^2);

% Integration limits
a = 0;
b = 2;

% Desired tolerance
tol = 1e-4;  % ~0.01% accuracy

% Initial step size
n = 2;  % start with 2 subintervals
h = (b - a)/n;

% Initial trapezoidal estimate
x = linspace(a, b, n+1);
S = sum(f(x)) - 0.5*(f(a)+f(b));
I_old = h*S;

iteration = 0;

while true
    iteration = iteration + 1;
    
    % Double the number of subintervals
    n = 2*n;
    h = (b - a)/n;
    x_new = a+h:h:b-h;   % new points only
    S_new = sum(f(x_new));
    
    % Trapezoidal estimate
    I_new = 0.5*I_old + h*S_new;  
    
    % Check tolerance
    if abs(I_new - I_old) <= tol
        break;
    end
    
    I_old = I_new;
end

fprintf('Adaptive Trapezoidal Integration:\n');
fprintf('Integral ≈ %.6f\n', I_new);
fprintf('Subintervals used: %d\n', n);
fprintf('Iterations: %d\n', iteration);
