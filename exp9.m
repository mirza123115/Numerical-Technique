% Euler Method for IVP: y' = (x - y)/2 , y(0) = 1
% Exact solution: y(x) = x - 2 + 3*exp(-x/2)

clear; clc; close all;

f = @(x,y) (x - y)/2;                % derivative function
y_exact = @(x) x - 2 + 3*exp(-x/2);  % exact solution

x0 = 0; 
y0 = 1; 
Xend = 3;

% Step sizes to test
h_values = [1, 0.5, 0.25, 0.125];
colors = lines(length(h_values));

figure; hold on; grid on;

% Plot exact solution (dense points for smooth curve)
x_exact = linspace(x0,Xend,500);
plot(x_exact, y_exact(x_exact),'k','LineWidth',2,'DisplayName','Exact');

% Loop over step sizes
for k = 1:length(h_values)
    h = h_values(k);
    N = (Xend-x0)/h;          
    x = zeros(1,N+1);
    y = zeros(1,N+1);
    
    % Initial values
    x(1) = x0; 
    y(1) = y0;
    
    % Euler method iteration
    for n = 1:N
        y(n+1) = y(n) + h*f(x(n),y(n));
        x(n+1) = x(n) + h;
    end
    
    % Plot numerical solution (smooth line only, no markers)
    plot(x,y,'-','Color',colors(k,:),'LineWidth',1.5,...
        'DisplayName',['h = ' num2str(h)]);
end

xlabel('x'); ylabel('y(x)');
title('Euler Method vs Exact Solution');
legend show; grid on;




%% RC charging: analytic + Euler numeric comparison


clear; clc; close all;

% Parameters
E = 20; R = 10; C = 117e-6; tau = R*C;
t_end = 3;

% Analytic solution (very fine grid)
t_fine = linspace(0, t_end, 20000);        % many points -> very smooth
Q_inf = C*E;
Q_exact = Q_inf * (1 - exp(-t_fine / tau));
i_exact = (E/R) * exp(-t_fine / tau);

% Euler numerical (fine step)
h = 1e-5;                                  % very small step -> smooth numeric
t_num = 0:h:t_end;
Q_num = zeros(size(t_num));
for k = 1:length(t_num)-1
    dQdt = (E - Q_num(k)/C)/R;
    Q_num(k+1) = Q_num(k) + h*dQdt;
end
i_num = (E - Q_num/C)/R;

% Plot
figure('Color','w','Units','normalized','Position',[0.1 0.1 0.6 0.6]);
subplot(2,1,1);
plot(t_fine, Q_exact,'k-','LineWidth',1.6); hold on;
plot(t_num, Q_num,'b-','LineWidth',1);      % no markers
xlabel('Time (s)'); ylabel('Q(t) [C]');
title('Capacitor Charge Q(t) - Analytic (black) & Euler (blue)');
xlim([0 0.01]);                             % zoom to transient
grid on;

subplot(2,1,2);
plot(t_fine, i_exact,'r-','LineWidth',1.6); hold on;
plot(t_num, i_num,'m-','LineWidth',1);
xlabel('Time (s)'); ylabel('i(t) [A]');
title('Current i(t) - Analytic (red) & Euler (magenta)');
xlim([0 0.01]);
grid on;






%% Improved Euler (Heun) vs Euler and exact solutions
clear; clc; close all;

%% ---------- Exercise 1: y' = (x-y)/2, y(0)=1 ----------
f1 = @(x,y) (x - y)/2;
y_exact = @(x) x - 2 + 3*exp(-x/2);    % given exact solution

x0 = 0; y0 = 1; Xend = 3;
h_values = [1, 0.5, 0.25, 0.125];

% figure with 2 subplots (solutions + error)
figure('Name','Exercise 1: Euler vs Improved Euler','Units','normalized','Position',[0.1 0.4 0.8 0.5]);
tiledlayout(1,2);

% --- subplot 1: curves ---
nexttile; hold on; grid on;
x_fine = linspace(x0,Xend,1000);
plot(x_fine, y_exact(x_fine),'k-','LineWidth',2,'DisplayName','Exact');

colors = lines(length(h_values));
results1 = [];

for k = 1:length(h_values)
    h = h_values(k);
    N = round((Xend - x0)/h);
    x = x0 + (0:N)*h;
    
    % Forward Euler
    y_eu = zeros(1,N+1); y_eu(1)=y0;
    for n=1:N
        y_eu(n+1) = y_eu(n) + h * f1(x(n), y_eu(n));
    end
    
    % Improved Euler (Heun)
    y_heun = zeros(1,N+1); y_heun(1)=y0;
    for n=1:N
        p = y_heun(n) + h * f1(x(n), y_heun(n)); % predictor
        y_heun(n+1) = y_heun(n) + (h/2) * ( f1(x(n), y_heun(n)) + f1(x(n+1), p) );
    end
    
    % plot curves
    plot(x, y_eu, '--','Color',colors(k,:), 'LineWidth',1.2, 'DisplayName',['Euler h=' num2str(h)]);
    plot(x, y_heun, '-','Color',colors(k,:), 'LineWidth',1.8, 'DisplayName',['Heun h=' num2str(h)]);
    
    % errors
    y_true_grid = y_exact(x);
    max_err_eu = max(abs(y_eu - y_true_grid));
    max_err_heun = max(abs(y_heun - y_true_grid));
    err_eu_3 = abs(y_eu(end) - y_exact(Xend));
    err_heun_3 = abs(y_heun(end) - y_exact(Xend));
    
    results1 = [results1; h, N, max_err_eu, max_err_heun, err_eu_3, err_heun_3]; %#ok<AGROW>
end

xlabel('x'); ylabel('y(x)');
title('Euler (dashed) & Heun (solid) vs Exact (black)');
legend('Location','eastoutside');

% --- subplot 2: error at x=3 ---
nexttile; hold on; grid on;
loglog(results1(:,1), results1(:,5), 'ro-','LineWidth',1.5,'MarkerFaceColor','r','DisplayName','Euler error');
loglog(results1(:,1), results1(:,6), 'bs-','LineWidth',1.5,'MarkerFaceColor','b','DisplayName','Heun error');
xlabel('Step size h'); ylabel('Error at x=3');
title('Convergence (error vs h)');
legend('Location','southwest');

% print error table
fprintf('\nExercise 1 errors (h, Nsteps, maxErr_Euler, maxErr_Heun, errAt3_Euler, errAt3_Heun):\n');
disp(results1);



%% ---------- Exercise 2: RC circuit ----------
% ODE: dQ/dt = (E - Q/C)/R
E = 20; R = 10; C = 117e-6;
tau = R*C; Q_inf = C*E;
t0 = 0; Q0 = 0; T_end = 3;

% analytic Q(t) and i(t)
Q_analytic = @(t) Q_inf * (1 - exp(-t/tau));
i_analytic = @(t) (E/R) * exp(-t/tau);

% choose time step for Euler/Heun for smooth numeric (small but not too small)
h_num = 1e-4;
Nnum = round((T_end - t0)/h_num);
t_num = t0 + (0:Nnum)*h_num;

% Forward Euler numeric
Q_eu = zeros(1,Nnum+1); Q_eu(1)=Q0;
for n=1:Nnum
    dQdt = (E - Q_eu(n)/C)/R;
    Q_eu(n+1) = Q_eu(n) + h_num * dQdt;
end
i_eu = (E - Q_eu./C) / R;

% Improved Euler (Heun) numeric
Q_heun = zeros(1,Nnum+1); Q_heun(1)=Q0;
for n=1:Nnum
    t_n = t_num(n);
    % predictor
    p = Q_heun(n) + h_num * ((E - Q_heun(n)/C)/R);
    % corrector
    f_n  = (E - Q_heun(n)/C)/R;
    f_np = (E - p/C)/R;
    Q_heun(n+1) = Q_heun(n) + (h_num/2) * (f_n + f_np);
end
i_heun = (E - Q_heun./C) / R;

% analytic arrays
t_fine = linspace(0, T_end, 20000);
Q_exact_vals = Q_analytic(t_fine);
i_exact_vals = i_analytic(t_fine);

% Plot RC results (zoom into transient region to see shape)
figure('Name','Exercise 2: RC circuit','Units','normalized','Position',[0.12 0.08 0.75 0.65]);
subplot(2,1,1);
plot(t_fine, Q_exact_vals, 'k-','LineWidth',1.6); hold on;
plot(t_num, Q_eu, 'r--','LineWidth',1);
plot(t_num, Q_heun, 'b-','LineWidth',1.2);
xlim([0, 0.01]); % transient region
xlabel('t (s)'); ylabel('Q(t) [C]'); title('Capacitor charge Q(t) (zoom 0..0.01 s)');
legend('Analytic', 'Euler','Heun','Location','best');

subplot(2,1,2);
plot(t_fine, i_exact_vals, 'k-','LineWidth',1.6); hold on;
plot(t_num, i_eu, 'r--','LineWidth',1);
plot(t_num, i_heun, 'b-','LineWidth',1.2);
xlim([0, 0.01]);
xlabel('t (s)'); ylabel('i(t) [A]'); title('Current i(t) (zoom 0..0.01 s)');
legend('Analytic', 'Euler','Heun','Location','best');

% Print a few numeric comparisons
fprintf('\nExercise 2: RC circuit (tau = %.3e s, Q_inf = %.6e C)\n', tau, Q_inf);
fprintf('At t = 1 ms, analytic Q = %.6e C, Euler Q = %.6e, Heun Q = %.6e\n', ...
    Q_analytic(1e-3), interp1(t_num,Q_eu,1e-3), interp1(t_num,Q_heun,1e-3));
fprintf('At t = 3 s, analytic Q = %.6e C, Euler Q = %.6e, Heun Q = %.6e\n', ...
    Q_analytic(3), Q_eu(end), Q_heun(end));

% End of script
