%% ================================================================
%  Creep nel modello di Maxwell–Zener (quasi-statico)
%  Deformazione ε(t) sotto sforzo costante σ0 H(t)
% ================================================================
clear; close all; clc;

%% Parametri
k0     = 10.0;
k      = 5.0;
beta   = 2.0;
sigma0 = 1.0;

tau_c = beta * (k0 + k) / (k0 * k);
fprintf('tau_c = %.4f\n', tau_c);

t_max = 6;          % tempo max (puoi adattarlo)
Nt    = 1000;
t     = linspace(0, t_max, Nt);

% Curva di creep
eps_t = (sigma0 / k0) * ( 1 - (k/(k0 + k)) .* exp(-t / tau_c) );
eps_0plus = sigma0 / (k0 + k);
eps_inf   = sigma0 / k0;

cols = lines(3);

fig = figure('Color','w'); hold on; grid on;
title('Curva di creep nel modello di Maxwell--Zener', ...
      'Interpreter','latex');
xlabel('$t$','Interpreter','latex','FontSize',14);
ylabel('$\varepsilon(t)$','Interpreter','latex','FontSize',14);

% Curva principale ε(t)
h_eps = plot(t, eps_t, 'LineWidth', 1.8, ...
             'Color', cols(1,:));

% Linea orizzontale per ε(∞) (arancione tratteggiata)
h_inf = yline(eps_inf, '--', ...
              'LineWidth', 1.4, ...
              'Color', cols(2,:), ...
              'HandleVisibility','off');

% Linea orizzontale per ε(0+) (gialla puntinata)
h_0   = yline(eps_0plus, ':', ...
              'LineWidth', 1.4, ...
              'Color', cols(3,:), ...
              'HandleVisibility','off');

% Marker iniziale
plot(0, eps_0plus, 'o', ...
     'MarkerSize', 6, ...
     'MarkerFaceColor', cols(3,:), ...
     'MarkerEdgeColor', 'k');

% Marker quasi asintotico (t ≈ 5 τ_c, se dentro l'intervallo)
t_mark = min(5 * tau_c, t_max);
[~, idx_mark] = min(abs(t - t_mark));
%plot(t(idx_mark), eps_t(idx_mark), 's', ...
%     'MarkerSize', 6, ...
%     'MarkerFaceColor', cols(2,:), ...
%     'MarkerEdgeColor', 'k');

% --- Limiti asse y (un po' sopra 0.1) ---
y_min = eps_0plus - 0.005;
y_max = eps_inf   + 0.005;
ylim([y_min, y_max]);
xlim([-0.2,6])

% --- Etichette testuali direttamente sul grafico ---
% posizione orizzontale delle scritte (fra 60% e 70% del grafico)
x_text_low  = 0.6 * t_max;
x_text_high = 0.6 * t_max;

text(x_text_low-2,  eps_0plus - 0.001, ...
     '$\frac{\sigma_0}{k_0+k}$', ...
     'Interpreter','latex', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','top', ...
     'Color','k');

text(x_text_high-2, eps_inf + 0.001, ...
     '$\frac{\sigma_0}{k_0}$', ...
     'Interpreter','latex', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'Color','k');

% Legenda solo per ε(t)
legend(h_eps, '$\varepsilon(t)$', ...
       'Interpreter','latex','Location','southeast');
