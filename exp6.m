%% Forward Difference Formula (O(h))

clc; clear;

f = @(x) exp(x);
x0 = 1;
exact = exp(1);

fprintf('h\t\tForward Diff\tError\n');
for k = 1:10
    h = 10^-k;
    f_prime = (f(x0+h) - f(x0))/h;
    error = abs(f_prime - exact);
    fprintf('%e\t%e\t%e\n', h, f_prime, error);
end


%% Central Difference Formula (O(h²))

clc; clear;

f = @(x) exp(x);
x0 = 1;
exact = exp(1);

fprintf('h\t\tCentral Diff O(h^2)\tError\n');
for k = 1:10
    h = 10^-k;
    f_prime = (f(x0+h) - f(x0-h))/(2*h);
    error = abs(f_prime - exact);
    fprintf('%e\t%e\t%e\n', h, f_prime, error);
end



%% Higher-Order Central Difference (O(h⁴))

clc; clear;

% Function
f = @(x) sin(cos(1./x));
x0 = 1;

% Initial step size and tolerance
h = 1;
tolerance = 1e-6;

% Variables to store previous results
D_prev2 = NaN;
D_prev = NaN;
iteration = 0;

fprintf('Iter\t h\t\tDerivative\tChange\n');

while true
    iteration = iteration + 1;
    
    % O(h^4) central difference formula
    D = (-f(x0 + 2*h) + 8*f(x0 + h) - 8*f(x0 - h) + f(x0 - 2*h)) / (12*h);
    
    % Compute change from previous iteration
    if ~isnan(D_prev)
        change = abs(D - D_prev);
    else
        change = NaN;
    end
    
    fprintf('%d\t%e\t%e\t%e\n', iteration, h, D, change);
    
    % Check convergence (stop if error stops decreasing or below tolerance)
    if iteration > 2
        if (change >= abs(D_prev - D_prev2)) || (change < tolerance)
            break;
        end
    end
    
    % Update previous results
    D_prev2 = D_prev;
    D_prev = D;
    
    % Reduce step size
    h = h / 10;
end

fprintf('\nApproximate derivative at x = %.2f is %.8f\n', x0, D);




%% Richardson's Extrapolation
clc; clear;

% Function
f = @(x) sin(x^3 - 7*x^2 + 6*x + 8);

x0 = 1.5;          % Point where derivative is needed
delta = 1e-13;     % Absolute error tolerance
toler = 1e-13;     % Relative error tolerance
maxIter = 12;      % Maximum iterations

% Initialization
err = 1; 
relerr = 1;
h = 1;
j = 1;
D = zeros(maxIter, maxIter);  % Richardson table

% First approximation using central difference
D(1,1) = (f(x0 + h) - f(x0 - h)) / (2*h);

fprintf('Iter\t h\t\tD_n\t\t\tError\t\tRelError\n');

while (relerr > toler) && (err > delta) && (j < maxIter)
    h = h / 2;
    % Central difference with smaller step
    D(j+1,1) = (f(x0 + h) - f(x0 - h)) / (2*h);
    
    % Richardson extrapolation
    for k = 1:j
        D(j+1,k+1) = D(j+1,k) + (D(j+1,k) - D(j,k)) / (4^k - 1);
    end
    
    % Compute error and relative error
    err = abs(D(j+1,j+1) - D(j,j));
    relerr = 2*err / (abs(D(j+1,j+1)) + abs(D(j,j)) + eps);
    
    fprintf('%d\t%e\t%e\t%e\t%e\n', j, h, D(j+1,j+1), err, relerr);
    
    j = j + 1;
end

fprintf('\nRichardson Approximation of f''(%.2f) is:\n%.15f\n', x0, D(j,j));
