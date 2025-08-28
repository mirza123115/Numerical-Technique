%% ===================== Exercise 1: Taylor Method 4th Order =====================
clear; clc; format long

a = 0; b = 3; y0 = 1; 
h_values = [1, 1/2, 1/4, 1/8];

f = @(x,y) x^2 - y; % ODE: y' = x^2 - y

fprintf('===== Exercise 1: Taylor 4th Order =====\n');
for idx = 1:length(h_values)
    h = h_values(idx);
    M = ceil((b-a)/h);
    T = a:h:b; T = T(:); % column vector
    Y_Taylor = zeros(M+1,1);
    Y_Taylor(1) = y0;
    
    for j = 1:M
        x = T(j); y = Y_Taylor(j);
        % Derivatives
        D1 = f(x,y);             % y'
        D2 = 2*x - D1;           % y''
        D3 = 2 - D2;             % y'''
        D4 = -D3;                % y''''
        % Taylor update
        Y_Taylor(j+1) = y + h*D1 + (h^2/2)*D2 + (h^3/6)*D3 + (h^4/24)*D4;
    end
    
    fprintf('Step size h = %g\n', h);
    disp(table(T, Y_Taylor, 'VariableNames', {'x','y_Taylor'}));
    figure; plot(T, Y_Taylor,'o-'); title(['Taylor Method, h = ', num2str(h)]); grid on;
end



%% Exercise 2: Taylor vs RK4 vs Exact Solution %%%%%%%%%%%%

clear; clc; format long;

a = 0; b = 1.4; ya = 0;
h = 0.01;  % small step
M = ceil((b-a)/h);
T = linspace(a,b,M+1)';

Y_RK = zeros(M+1,1);
Y_RK(1) = ya;

f_fun = @(t,y) t^2 - y;

for j = 1:M
    k1 = h*f_fun(T(j), Y_RK(j));
    k2 = h*f_fun(T(j)+h/2, Y_RK(j)+k1/2);
    k3 = h*f_fun(T(j)+h/2, Y_RK(j)+k2/2);
    k4 = h*f_fun(T(j)+h, Y_RK(j)+k3);
    Y_RK(j+1) = Y_RK(j) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

Y_exact = exp(-1) + T.^2 - 2*T + 2;

% Plot
figure;
plot(T,Y_exact,'r--','LineWidth',1.5); hold on;
plot(T,Y_RK,'b-','LineWidth',1.5);
xlabel('t'); ylabel('y(t)');
legend('Exact','RK4');
grid on;








%% Exercise 3: RK4 for 2nd-Order ODE %%%%%%%%%
clear; clc; format long

% Endpoints and initial conditions
a = 0; b = 5;
Za = [3 -5]; % x(0)=3, x'(0)=-5

% Step size and number of steps
M = 50;
h = (b-a)/M;

% Time vector
T = a:h:b; 
T = T(:); % column vector

% Initialize solution matrix: Z(:,1)=x, Z(:,2)=x')
Z = zeros(M+1, length(Za));
Z(1,:) = Za;

% RK4 method
for j=1:M
    k1 = h*F3(T(j), Z(j,:));
    k2 = h*F3(T(j)+h/2, Z(j,:)+k1/2);
    k3 = h*F3(T(j)+h/2, Z(j,:)+k2/2);
    k4 = h*F3(T(j)+h, Z(j,:)+k3);
    
    Z(j+1,:) = Z(j,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

% Numerical solution
x_num = Z(:,1);

% Exact solution
x_exact = exp(-2*T).*(3*cos(T) + sin(T));

% Display results in table
disp(table(T, x_num, x_exact, 'VariableNames', {'t','x_RK4','x_exact'}));

% Plot numerical vs exact solution
figure;
plot(T, x_num, 'bo-', T, x_exact, 'r--','LineWidth',1.5);
xlabel('t'); ylabel('x(t)');
legend('RK4 Numerical','Exact Solution');
title('Exercise 3: RK4 vs Exact Solution for x'''' + 4x'' + 5x = 0');
grid on;

%%%%%%%%% Function for the system of ODEs %%%%%%%%%
% function Zout = F3(t,Z)
%     x = Z(1);
%     y = Z(2);
%     Zout = [y, -5*x - 4*y]; % [x', y']
% end



%% Exercise 4: RK4 Linear Shooting Method + Exact Solution %%%%%%%%
clear; clc; format long

% Interval and steps
a = 0; b = 4; h = 0.1; 
M = ceil((b-a)/h); 
T = a:h:b; 
T = T(:);  % column vector

% Boundary values
alpha = 1.25; 
beta  = -0.95;

% Initialize solution vectors
u = zeros(M+1,1); up = zeros(M+1,1);
v = zeros(M+1,1); vp = zeros(M+1,1);

% Initial conditions for the two IVPs
u(1) = alpha; up(1) = 0;
v(1) = 0;    vp(1) = 1;

% RHS function f(t,x,x')
F = @(t,x,xd) (2*t/(1+t^2))*xd - (2/(1+t^2))*x + 1;

% Solve u(t) with RK4
for j = 1:M
    k1 = up(j); l1 = F(T(j), u(j), up(j));
    k2 = up(j)+h*l1/2; l2 = F(T(j)+h/2, u(j)+h*k1/2, up(j)+h*l1/2);
    k3 = up(j)+h*l2/2; l3 = F(T(j)+h/2, u(j)+h*k2/2, up(j)+h*l2/2);
    k4 = up(j)+h*l3;   l4 = F(T(j)+h, u(j)+h*k3, up(j)+h*l3);
    
    u(j+1)  = u(j)  + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    up(j+1) = up(j) + h*(l1 + 2*l2 + 2*l3 + l4)/6;
end

% Solve v(t) with RK4
for j = 1:M
    k1 = vp(j); l1 = F(T(j), v(j), vp(j));
    k2 = vp(j)+h*l1/2; l2 = F(T(j)+h/2, v(j)+h*k1/2, vp(j)+h*l1/2);
    k3 = vp(j)+h*l2/2; l3 = F(T(j)+h/2, v(j)+h*k2/2, vp(j)+h*l2/2);
    k4 = vp(j)+h*l3;   l4 = F(T(j)+h, v(j)+h*k3, vp(j)+h*l3);
    
    v(j+1)  = v(j)  + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    vp(j+1) = vp(j) + h*(l1 + 2*l2 + 2*l3 + l4)/6;
end

% Compute C for linear combination
C = (beta - u(end)) / v(end);

% Final numerical solution of BVP
x_bvp = u + C*v;

% Exact solution
x_exact = 1.25 + 0.4860896526*T - 2.25*T.^2 + 2*T.*atan(T) ...
          - 0.5*log(1+T.^2) + 0.5*T.^2.*log(1+T.^2);

% Compute error
error = abs(x_bvp - x_exact);

% Display results in a table
disp(table(T, x_bvp, x_exact, error, 'VariableNames', {'t','x_BVP','x_exact','Error'}));

% Plot numerical vs exact solution
figure;
plot(T, x_bvp, 'b-o','LineWidth',1.5); hold on;
plot(T, x_exact, 'r--','LineWidth',1.5);
xlabel('t'); ylabel('x(t)');
legend('RK4 BVP','Exact Solution');
title('Exercise 4: RK4 Linear Shooting BVP vs Exact Solution');
grid on;

% Plot absolute error
figure;
plot(T, error, 'k','LineWidth',1.5);
xlabel('t'); ylabel('Absolute Error');
title('Absolute Error between RK4 Solution and Exact Solution');
grid on;
