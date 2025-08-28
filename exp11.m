clear all
close all
clc
%% Exercise 1: Laplace Equation (Liebmann's method)
% Plate: 4x4 cm, dx = dy = 1 cm
% Example boundary: top=1000, bottom=0, left=750, right=500 (adjust if needed)



% Grid size
nx = 5; ny = 5; % 0:4 in both directions -> 5 nodes each
T = zeros(ny, nx);

% Boundary conditions
T(1,:)   = 0;     % bottom edge y=0
T(end,:) = 1000;  % top edge y=4
T(:,1)   = 750;   % left edge x=0
T(:,end) = 500;   % right edge x=4

% Iteration parameters
tol = 0.01; % 1% error tolerance
error = 1; 
iter = 0;

while error > tol
    error = 0;
    for j = 2:ny-1
        for i = 2:nx-1
            old = T(j,i);
            T(j,i) = 0.25*(T(j+1,i) + T(j-1,i) + T(j,i+1) + T(j,i-1));
            error = max(error, abs((T(j,i)-old)/T(j,i))*100);
        end
    end
    iter = iter + 1;
end

fprintf('Converged in %d iterations\n', iter);
disp(T);

% Plot
[X,Y] = meshgrid(0:4,0:4);
figure; surf(X,Y,T);
xlabel('x (cm)'); ylabel('y (cm)'); zlabel('Temperature (°C)');
title('2D Laplace Equation (Liebmann Method)');





%% Exercise 2: 1D Heat equation explicit method (FTCS)
clear; clc;

L = 10; dx = 2; x = 0:dx:L; 
dt = 1.0; 
tmax = 4;
r = dt/dx^2; % with k=1 normalized

% Initial condition
T = zeros(size(x));
T(1) = 1000; % left boundary
T(end) = 500; % right boundary

nsteps = tmax/dt;

for n = 1:nsteps
    Tn = T;
    for j = 2:length(x)-1
        T(j) = Tn(j) + r*(Tn(j+1)-2*Tn(j)+Tn(j-1));
    end
end

disp('Explicit solution at t=4 s:');
disp(table(x',T','VariableNames',{'x_cm','Temperature'}));

plot(x,T,'-o'); grid on;
xlabel('x (cm)'); ylabel('Temperature (°C)');
title('Explicit FTCS Heat Equation at t=4 s');




%% Exercise 3: 1D Heat equation implicit method (Backward Euler)

clear; clc;

L = 10; dx = 2; x = 0:dx:L; 
dt = 1.0; 
tmax = 4;
r = dt/dx^2; 

% Initial condition
T = zeros(length(x),1);  % make column vector
T(1) = 1000; 
T(end) = 500;

nsteps = tmax/dt;
N = length(x);

% Construct coefficient matrix A (tridiagonal)
A = (1+2*r)*eye(N);
for i=2:N-1
    A(i,i-1) = -r;
    A(i,i+1) = -r;
end
% Boundary conditions: keep T(1), T(N) fixed
A(1,:) = 0; A(1,1)=1;
A(N,:) = 0; A(N,N)=1;

for n=1:nsteps
    % Right-hand side is old T
    b = T;  
    % Enforce boundary values
    b(1) = 1000; 
    b(N) = 500;
    % Solve system
    T = A\b;
end

disp('Implicit solution at t=4 s:');
disp(table(x',T,'VariableNames',{'x_cm','Temperature'}));

plot(x,T,'-o'); grid on;
xlabel('x (cm)'); ylabel('Temperature (°C)');
title('Implicit Backward Euler Heat Equation at t=4 s');
