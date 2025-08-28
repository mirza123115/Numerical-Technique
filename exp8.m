clc; clear;

% Define the function
f = @(x) 5*x + 1 + x^3;

% Bisection method parameters
xLow = -1;
xUp = 0;
xTol = 4e-4;

% Initial midpoint
xMid = (xLow + xUp)/2;
yMid = f(xMid);

iters = 0;

fprintf('Iter\t xLow\t xUp\t xMid\t f(xMid)\n');

while ( (xUp - xLow)/2 > xTol )
    iters = iters + 1;
    
    yLow = f(xLow);
    
    % Check which subinterval contains the root
    if yLow * yMid > 0
        xLow = xMid;
    else
        xUp = xMid;
    end
    
    % Compute new midpoint
    xMid = (xLow + xUp)/2;
    yMid = f(xMid);
    
    fprintf('%d\t%e\t%e\t%e\t%e\n', iters, xLow, xUp, xMid, yMid);
end

fprintf('\nApproximate root: %f\n', xMid);
fprintf('Function value at root: %e\n', yMid);
fprintf('Iterations: %d\n', iters);



%% False-Position Method (Regula Falsi)

clc; clear;

f = @(x) 5*x^3 + x - 1; % example function
xLow = -1;
xUp = 0;
epsilon = 1e-4;

iteration = 0;

while true
    iteration = iteration + 1;
    
    fLow = f(xLow);
    fUp = f(xUp);
    
    % False position formula
    xNew = xUp - fUp*(xUp - xLow)/(fUp - fLow);
    fNew = f(xNew);
    
    fprintf('Iter %d: x = %.6f, f(x) = %.6f\n', iteration, xNew, fNew);
    
    % Check convergence
    if abs(fNew) < epsilon
        break;
    end
    
    % Update bracket
    if fLow * fNew < 0
        xUp = xNew;
    else
        xLow = xNew;
    end
end

fprintf('Approximate root: %.6f after %d iterations\n', xNew, iteration);


%% Newton-Raphson Method

clc; clear;

% Function and derivative
f = @(x) 1 - x.*exp(x);
fDeriv = @(x) -(1 + x).*exp(x);

% Initial guess and tolerance
x = 0;          
tol = 8e-8;     
iteration = 0;

fprintf('Iter\t x_k\t\tf(x_k)\n');

while true
    iteration = iteration + 1;
    
    dx = -f(x)/fDeriv(x);   % Newton-Raphson formula
    x = x + dx;
    
    fprintf('%d\t %.10f\t %.10f\n', iteration, x, f(x));
    
    if abs(dx) < tol
        break;
    end
end

fprintf('Estimated root: %.10f after %d iterations\n', x, iteration);


%% Secant Method


clc; clear;

% Function definition
f = @(x) x.^3 + sin(x) - exp(x);

% Initial guesses
x0 = 0;
x1 = 1;

% Tolerance and maximum iterations
tol = 1e-7;
maxIter = 100;

% Initial function values
f0 = f(x0);
f1 = f(x1);

iteration = 0;

fprintf('Iter\t x_k\t\tf(x_k)\n');

while abs(x1 - x0) > tol && iteration < maxIter
    iteration = iteration + 1;
    
    % Secant formula
    x_new = x1 - f1*(x1 - x0)/(f1 - f0);
    
    fprintf('%d\t %.10f\t %.10f\n', iteration, x_new, f(x_new));
    
    % Update for next iteration
    x0 = x1;
    f0 = f1;
    
    x1 = x_new;
    f1 = f(x1);
end

fprintf('Estimated root: %.10f after %d iterations\n', x1, iteration);
