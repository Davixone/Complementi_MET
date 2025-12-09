%% ================================================================
%  Studio di orbite nello spazio delle fasi
%  ================================================================

clear; close all; clc;

%% ---------------------------- PARAMETRI FISICI ----------------------------
m    = 1.0;      % massa
k0   = 10.0;     % molla esterna
k    = 5.0;      % molla del ramo di Maxwell
beta = 2.0;      % viscosità del dashpot

params = struct('m',m,'k0',k0,'k',k,'beta',beta);

%% --------------------- INTERVALLO E CONDIZIONI INIZIALI --------------------
tspan = [0 20];

ICs = [
    1   0   0;     % x(0)=1, xdot(0)=0, x1(0)=0
    0   1   0;     % x(0)=0, xdot(0)=1
    1  -1   0;     % x(0)=1, xdot(0)=-1
    0.5 0  0.5     % x(0)=0.5, x1(0)=0.5
];

colors = lines(size(ICs,1));

%% ========================================================================
%  FIGURA 1 — Orbite 3D nello spazio delle fasi (x, xdot, x1)
% ========================================================================
fig1 = figure('Color','w');  % <-- sfondo bianco
hold on; grid on;

title('Orbite nello spazio delle fasi $(x,\dot{x},x_1)$','Interpreter','latex');
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',14);
zlabel('$x_1$','Interpreter','latex','FontSize',14);

view(135,30);

for i = 1:size(ICs,1)
    z0 = ICs(i,:).';
    ode = @(t,z) maxwell_system(t,z,params);
    [t,z] = ode45(ode, tspan, z0);

    plot3(z(:,1), z(:,2), z(:,3), ...
          'Color', colors(i,:), 'LineWidth', 1.5);
end

plot3(0,0,0,'ko','MarkerFaceColor','k'); % equilibrio

legend('IC1','IC2','IC3','IC4','Equilibrio','Location','best','Interpreter','latex');

%% ========================================================================
%  FIGURA 2 — Proiezione 2D nello spazio delle fasi (x, xdot)
% ========================================================================
fig2 = figure('Color','w');
hold on; grid on;

title('Proiezione nello spazio delle fasi $(x,\dot{x})$','Interpreter','latex');
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',14);

for i = 1:size(ICs,1)
    z0 = ICs(i,:).';
    ode = @(t,z) maxwell_system(t,z,params);
    [t,z] = ode45(ode, tspan, z0);

    plot(z(:,1), z(:,2), 'Color', colors(i,:), 'LineWidth', 1.5);
end

plot(0,0,'ko','MarkerFaceColor','k');
legend('IC1','IC2','IC3','IC4','Equilibrio','Location','best','Interpreter','latex');

%% ========================================================================
%  FIGURA 3 — Effetto della viscosità beta sulle orbite
% ========================================================================
betas = [0.4 2 10];   % viscosità bassa, media, alta

% Colori originali 'lines'
colors2 = lines(numel(betas));

% Stili delle linee:
% blu piena, arancione tratteggiata, gialla piena
lineStyles = {'-', '--', ':'};

fig3 = figure('Color','w');
hold on; grid on;

title('Effetto di $\beta$ sulle orbite $(x,\dot{x})$','Interpreter','latex');
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',14);

for i = 1:numel(betas)
    params.beta = betas(i);
    z0 = [1; 0; 0]; 
    ode = @(t,z) maxwell_system(t,z,params);
    [t,z] = ode45(ode, tspan, z0);
    
    plot(z(:,1), z(:,2), ...
         'Color', colors2(i,:), ...
         'LineStyle', lineStyles{i}, ...
         'LineWidth', 1.8);
end
plot(0,0,'ko','MarkerFaceColor','k'); % Equilibrio al centro
legend('$\beta=0.4$','$\beta=2$','$\beta=10$', 'Equilibrio', ...
       'Interpreter','latex','Location','best');


%% ========================================================================
%  FIGURA 4 — Effetto della rigidezza k (ramo di Maxwell) sulle orbite
% ========================================================================
k_values = [4 10 50];       % k piccolo, intermedio, grande
colors_k = lines(numel(k_values));
lineStyles_k = {'-','--',':'};

fig4 = figure('Color','w');
hold on; grid on;

title('Effetto di $k$ sulle orbite $(x,\dot{x})$','Interpreter','latex');
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',14);

for i = 1:numel(k_values)
    params.k = k_values(i);       % varia solo k
    params.k0 = 10;               % k0 fissato
    params.beta = 2;              % beta fissato
    params.m = 1;                 % m fissata
    
    z0 = [1; 0; 0];               % stessa condizione iniziale
    ode = @(t,z) maxwell_system(t,z,params);
    [t,z] = ode45(ode, tspan, z0);

    plot(z(:,1), z(:,2), ...
         'Color', colors_k(i,:), ...
         'LineStyle', lineStyles_k{i}, ...
         'LineWidth', 1.8);
end

h_eq = plot(0,0,'ko','MarkerFaceColor','k'); % equilibrio
legend({'$k=4$','$k=10$','$k=50$','Equilibrio'}, ...
       'Interpreter','latex','Location','best');


%% ========================================================================
%  FIGURA 5 — Effetto della rigidezza k0 (molla esterna) sulle orbite
% ========================================================================
k0_values = [4 10 50];     % k0 piccolo, intermedio, grande
colors_k0 = lines(numel(k0_values));
lineStyles_k0 = {'-','--',':'};

fig5 = figure('Color','w');
hold on; grid on;

title('Effetto di $k_0$ sulle orbite $(x,\dot{x})$','Interpreter','latex');
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$\dot{x}$','Interpreter','latex','FontSize',14);

for i = 1:numel(k0_values)
    params.k0 = k0_values(i);     % varia solo k0
    params.k  = 5;                % k fissato
    params.beta = 2;              % beta fissato
    params.m = 1;                 % m fissata
    
    z0 = [1; 0; 0];               % stessa condizione iniziale
    ode = @(t,z) maxwell_system(t,z,params);
    [t,z] = ode45(ode, tspan, z0);

    plot(z(:,1), z(:,2), ...
         'Color', colors_k0(i,:), ...
         'LineStyle', lineStyles_k0{i}, ...
         'LineWidth', 1.8);
end

h_eq = plot(0,0,'ko','MarkerFaceColor','k'); % equilibrio
legend({'$k_0=4$','$k_0=10$','$k_0=50$','Equilibrio'}, ...
       'Interpreter','latex','Location','best');




%% ========================================================================
%  FUNZIONE — Sistema dinamico di Maxwell
% ========================================================================
function dz = maxwell_system(~, z, p)
    % Variabili di stato
    x   = z(1);
    xd  = z(2);
    x1  = z(3);

    % Parametri
    m    = p.m;
    k0   = p.k0;
    k    = p.k;
    beta = p.beta;

    % Sistema dinamico
    dz   = zeros(3,1);
    dz(1) = xd;
    dz(2) = - (k0/m)*x - (k/m)*x1;
    dz(3) = xd - (k/beta)*x1;
end


%% ================================================================
%  Confronto soluzione analitica vs numerica per x(t)
%  ================================================================
clear; close all; clc;

%% Parametri fisici
m    = 1.0;
k0   = 10.0;
k    = 5.0;
betas = [0.5 2 10];   % tre regimi di smorzamento

x0   = 1.0;           % condizioni iniziali
v0   = 0.0;
x10  = 0.0;
z0   = [x0; v0; x10];

tspan = [0 20];       % intervallo temporale per integrazione
t_plot = linspace(0,20,1000);   % tempi per valutare la soluzione analitica

colors = lines(numel(betas));
lineStyles = {'-','-','-'};     % linee continue per l'analitica
markers    = {'o','s','^'};     % marcatori per la numerica

figure('Color','w'); hold on; grid on;
title('Confronto tra soluzione analitica e numerica per $x(t)$', ...
      'Interpreter','latex');
xlabel('$t$','Interpreter','latex','FontSize',14);
ylabel('$x(t)$','Interpreter','latex','FontSize',14);

for i = 1:numel(betas)
    beta = betas(i);
    
    %----------------- Matrice A per questo beta -----------------
    A = [ 0         1          0;
         -k0/m      0        -k/m;
          0         1      -k/beta ];
    
    %----------------- Soluzione analitica via expm -----------------
    x_analitica = zeros(size(t_plot));
    for j = 1:numel(t_plot)
        z_t = expm(A * t_plot(j)) * z0;
        x_analitica(j) = z_t(1);
    end
    
    %----------------- Soluzione numerica con ode15s -----------------
    params = struct('m',m,'k0',k0,'k',k,'beta',beta);
    ode = @(t,z) maxwell_system(t,z,params);
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t_num,z_num] = ode15s(ode, tspan, z0, opts);
    
    %----------------- Plot: analitica (linea) vs numerica (punti) ---
    plot(t_plot, x_analitica, ...
         'Color', colors(i,:), ...
         'LineStyle', lineStyles{i}, ...
         'LineWidth', 1.6);
    
    plot(t_num, z_num(:,1), ...
         markers{i}, ...
         'Color', colors(i,:), ...
         'MarkerSize', 4, ...
         'MarkerFaceColor', colors(i,:), ...
         'LineStyle', 'none');
end

legend({'Analitica $\beta=0.5$','Numerica $\beta=0.5$', ...
        'Analitica $\beta=2$','Numerica $\beta=2$', ...
        'Analitica $\beta=10$','Numerica $\beta=10$'}, ...
        'Interpreter','latex','Location','best');

%% ================================================================
%  Oscillatore di Maxwell–Zener
%  x(t) analitica + inviluppi ∝ e^{-alpha t}
%  ================================================================
clear; close all; clc;

%% Parametri fisici
m    = 1.0;
k0   = 10.0;
k    = 5.0;
betas = [0.5 10];

% Condizioni iniziali
x0   = 1.0;
v0   = 0.0;
x10  = 0.0;
z0   = [x0; v0; x10];

% Tempo
t_max  = 20;
Nt     = 1000;
t_plot = linspace(0, t_max, Nt);

colors = lines(numel(betas));

figure('Color','w'); hold on; grid on;
title('Risposta $x(t)$ con inviluppi proporzionali a $e^{-\alpha t}$', ...
      'Interpreter','latex');
xlabel('$t$','Interpreter','latex','FontSize',14);
ylabel('$x(t)$','Interpreter','latex','FontSize',14);

legend_handles = [];
legend_labels  = {};

for i = 1:numel(betas)

    beta = betas(i);

    % Matrice dinamica
    A = [ 0         1          0;
         -k0/m      0        -k/m;
          0         1      -k/beta ];

    % Autovalori → alpha
    lambda = eig(A);
    lambda_c = lambda(abs(imag(lambda)) > 1e-8);
    alpha = -real(lambda_c(1));
    fprintf('beta = %.2f -> alpha = %.4f\n', beta, alpha);

    % Soluzione analitica via expm
    x_an = zeros(size(t_plot));
    for j = 1:numel(t_plot)
        z_t     = expm(A * t_plot(j)) * z0;
        x_an(j) = z_t(1);
    end

    % ======= traccia x(t) =======
    h_main = plot(t_plot, x_an, ...
         'Color', colors(i,:), 'LineWidth', 1.6);
    legend_handles(end+1) = h_main;
    legend_labels{end+1}  = sprintf('$x(t),\\ \\beta=%.1f$', beta);

    % ======= inviluppi =======
    E   = max(abs(x_an));
    env = E * exp(-alpha * t_plot);

    plot(t_plot,  env, '--', 'Color', colors(i,:), ...
         'LineWidth', 1.2, 'HandleVisibility','off');
    plot(t_plot, -env, '--', 'Color', colors(i,:), ...
         'LineWidth', 1.2, 'HandleVisibility','off');

    % Voce singola in legenda per il tratteggio (senza raddoppiare)
    h_env = plot(nan, nan, '--', 'Color', colors(i,:), 'LineWidth', 1.2);
    legend_handles(end+1) = h_env;
    legend_labels{end+1}  = sprintf('$\\pm E e^{-\\alpha t},\\ \\beta=%.1f$', beta);
end

legend(legend_handles, legend_labels, ...
    'Interpreter','latex', 'Location','best');
