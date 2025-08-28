% Coefficient matrix
A = [2 1 -1;
    -3 -1 2;
    -2 1 2];

% Right-hand side
b = [3; -3; -3];

% Using backslash operator
x1 = A\b;
disp('Solution using backslash operator:');
disp(x1)

% Using pseudo-inverse
x2 = pinv(A)*b;
disp('Solution using pseudo-inverse:');
disp(x2)





%% Gaussian Elimination with Back Substitution

clc; clear;

% Coefficient matrix A
A = [3 4 1;
     6 7 2;
     2 3 5];

% Right-hand side vector b
b = [14; 26; 23];

% Augmented matrix [A|b]
Aug = [A b];

[n, ~] = size(Aug);

%% Triangularization (Forward Elimination)
for k = 1:n-1
    for i = k+1:n
        factor = Aug(i,k)/Aug(k,k);
        Aug(i,k:end) = Aug(i,k:end) - factor*Aug(k,k:end);
    end
end

disp('Upper Triangular Form [A|b]:');
disp(Aug);

%% Back Substitution
x = zeros(n,1);
x(n) = Aug(n,end)/Aug(n,n);

for i = n-1:-1:1
    sum = 0;
    for j = i+1:n
        sum = sum + Aug(i,j)*x(j);
    end
    x(i) = (Aug(i,end) - sum)/Aug(i,i);
end

disp('Solution [x1; x2; x3]:');
disp(x);




%% MATLAB Code: Gaussian Elimination with Pivoting

clc; clear;

% Example system: modify A and b for any system
A = [2 4 -6;
     1 1 6;
     -1 2 9];

b = [-4; 7; 2];

[n, ~] = size(A);

Aug = [A b];

%% Gauss-Jordan Elimination with Pivoting
for k = 1:n
    % --- Pivoting ---
    [~, pk] = max(abs(Aug(k:n,k)));
    pk = pk + k - 1;
    if pk ~= k
        Aug([k pk], :) = Aug([pk k], :); % swap rows
    end
    
    % --- Make pivot 1 ---
    Aug(k,:) = Aug(k,:) / Aug(k,k);
    
    % --- Eliminate other entries in column k ---
    for i = 1:n
        if i ~= k
            factor = Aug(i,k);
            Aug(i,:) = Aug(i,:) - factor*Aug(k,:);
        end
    end
end

disp('Reduced Row Echelon Form [A|b]:');
disp(Aug);

solution = Aug(:,end);
disp('Solution [x1; x2; x3]:');
disp(solution);





%% Gauss-Seidel Iterative Method

clc; clear;

% Initial guess
x = [0;0;0];
tol = 1e-6;
maxIter = 100;

for k = 1:maxIter
    x_old = x;
    
    % Gauss-Seidel updates
    x(1) = (12 - x(2) - x(3))/4;
    x(2) = (20 - 2*x(1) - x(3))/7;
    x(3) = (30 - x(1) - 2*x(2))/8;
    
    %Check convergence
    if norm(x - x_old, inf) < tol
        break;
    end
end

disp(['Gauss-Seidel solution after ', num2str(k), ' iterations:']);
disp(x);
