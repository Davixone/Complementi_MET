%% ================================================================
%  Relaxation nel modello di Maxwell–Zener (quasi-statico)
%  Tensione σ(t) sotto deformazione imposta ε0 Θ(t)
% ================================================================
clear; close all; clc;

%% Parametri del modello
k0    = 10.0;   % molla elastica pura
k     = 5.0;    % molla nel ramo di Maxwell
beta  = 2.0;    % viscosità del dashpot
eps0  = 0.01;   % deformazione imposta (piccola, ad es. 1%)

% Tempo caratteristico di relaxation
tau_r = beta / k;
fprintf('tau_r = %.4f\n', tau_r);

%% Intervallo temporale
t_max = 6;          % puoi adattarlo
Nt    = 1000;
t     = linspace(0, t_max, Nt);

%% Curva di relaxation σ(t)
% σ(t) = k0*ε0 + k*ε0 * exp(-t/τ_r)
sigma_t   = k0*eps0 + k*eps0 .* exp(-t / tau_r);
sigma_0p  = (k0 + k) * eps0;   % σ(0+)
sigma_inf = k0 * eps0;         % σ(∞)

%% Colori e figura
cols = lines(3);

fig = figure('Color','w'); hold on; grid on;
title('Curva di relaxation nel modello di Maxwell--Zener', ...
      'Interpreter','latex');
xlabel('$t$','Interpreter','latex','FontSize',14);
ylabel('$\sigma(t)$','Interpreter','latex','FontSize',14);

% Curva principale σ(t)
h_sig = plot(t, sigma_t, 'LineWidth', 1.8, ...
             'Color', cols(1,:));

% Linea orizzontale per σ(∞) (arancione tratteggiata)
yline(sigma_inf, '--', 'LineWidth', 1.4, ...
      'Color', cols(2,:), 'HandleVisibility','off');

% Linea orizzontale per σ(0+) (gialla puntinata)
yline(sigma_0p, ':', 'LineWidth', 1.4, ...
      'Color', cols(3,:), 'HandleVisibility','off');

% Marker al tempo t=0+
plot(0, sigma_0p, 'o', ...
     'MarkerSize', 6, ...
     'MarkerFaceColor', cols(3,:), ...
     'MarkerEdgeColor', 'k');

% Marker su un punto "quasi rilassato", ad esempio t ≈ 5 τ_r
t_mark = min(5 * tau_r, t_max);
[~, idx_mark] = min(abs(t - t_mark));
%plot(t(idx_mark), sigma_t(idx_mark), 's', ...
%     'MarkerSize', 6, ...
%     'MarkerFaceColor', cols(2,:), ...
%     'MarkerEdgeColor', 'k');

% Limiti asse y (leggermente più larghi del range)
y_min = sigma_inf  - 0.1 * (sigma_0p - sigma_inf);
y_max = sigma_0p   + 0.1 * (sigma_0p - sigma_inf);
ylim([y_min, y_max]);
xlim([-0.1,6]);

% Etichette testuali direttamente sul grafico
x_text_low  = 0.6 * t_max;
x_text_high = 0.6 * t_max;

text(x_text_high, sigma_0p + 0.01*(sigma_0p - sigma_inf), ...
     '$(k_0 + k)\,\varepsilon_0$', ...
     'Interpreter','latex', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'Color','k');

text(x_text_low,  sigma_inf - 0.01*(sigma_0p - sigma_inf), ...
     '$k_0\,\varepsilon_0$', ...
     'Interpreter','latex', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','top', ...
     'Color','k');

% Legenda solo per σ(t)
legend(h_sig, '$\sigma(t)$', ...
       'Interpreter','latex','Location','northeast');

% Salvataggio figura
set(fig,'PaperPositionMode','auto');
print(fig,'relaxation_curve_maxwell_zener.png','-dpng','-r400');
