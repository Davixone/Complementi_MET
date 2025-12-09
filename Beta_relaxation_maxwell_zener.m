%% ================================================================
%  Relaxation nel modello di Maxwell–Zener
%  Variazione di β: tre curve (blu, arancione, giallo)
% ================================================================
clear; close all; clc;

%% Parametri base
k0 = 10;
k  = 5;
eps0 = 0.01;   % deformazione imposta

betas = [0.5 2 10];

colors = lines(numel(betas));
styles = {'-','--',':'};

%% Tempo
tmax = 10;
t = linspace(0,tmax,1000);

fig = figure('Color','w'); hold on; grid on;
title('Relaxation: variazione della tensione $\sigma(t)$ al variare di $\beta$', ...
      'Interpreter','latex');
xlabel('$t$','Interpreter','latex','FontSize',14);
ylabel('$\sigma(t)$','Interpreter','latex','FontSize',14);
xlim([-0.1,10]);
ylim([0.09,0.16]);



for i = 1:numel(betas)
    beta = betas(i);

    % tempo caratteristico di relaxation
    tau_r = beta / k;

    % soluzione analitica della relaxation
    sig_t = k0*eps0 + k*eps0 * exp(-t/tau_r);

    plot(t, sig_t, ...
        'Color', colors(i,:), ...
        'LineStyle', styles{i}, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('$\\beta = %.1f$', beta));
end

legend('Interpreter','latex','Location','best');
