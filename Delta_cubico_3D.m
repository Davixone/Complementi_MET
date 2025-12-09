    %% ================================================================
    %  Discriminante del cubico λ^3 + a2 λ^2 + a1 λ + a0
    %  in funzione di (a2, a1, a0) nel primo ottante
    %  ================================================================
    
    clear; close all; clc;
    
    % ----------------- Funzione discriminante per un cubico -----------
    % Per p(λ) = λ^3 + b λ^2 + c λ + d  (a=1) vale:
    % Δ = 18 b c d - 4 b^3 d + b^2 c^2 - 4 c^3 - 27 d^2
    cubic_discriminant = @(b,c,d) ...
        18.*b.*c.*d - 4.*(b.^3).*d + (b.^2).*(c.^2) ...
        - 4.*(c.^3) - 27.*(d.^2);
    
    % ----------------- Dominio dei coefficienti -----------------------
    a2_vec = linspace(0.1, 10, 20);   % coefficiente di λ^2
    a1_vec = linspace(0.1, 10, 20);   % coefficiente di λ
    a0_vec = linspace(0.1, 10, 20);   % termine noto
    
    [A2, A1, A0] = ndgrid(a2_vec, a1_vec, a0_vec);
    
    % ----------------- Discriminante su tutta la griglia --------------
    Delta = cubic_discriminant(A2, A1, A0);
    
    % ----------------- Sottocampionamento per lo scatter --------------
    Npoints = 2500;  % riduci/aumenta a piacere
    idx = randperm(numel(Delta), Npoints);
    
    x = A2(idx);   % asse x: a2
    y = A1(idx);   % asse y: a1
    z = A0(idx);   % asse z: a0
    dvals = Delta(idx);
    
    % ----------------- Grafico 3D con colore = Δ ----------------------
    fig = figure('Color','w');
    ax  = axes(fig); hold(ax,'on'); grid(ax,'on');
    sc  = scatter3(ax, x, y, z, 15, dvals, 'filled');
    
    xlabel('$a_2$','Interpreter','latex','FontSize',13);
    ylabel('$a_1$','Interpreter','latex','FontSize',13);
    zlabel('$a_0$','Interpreter','latex','FontSize',13);
    title('$\Delta$ del polinomio $\lambda^3 + a_2 \lambda^2 + a_1 \lambda + a_0$', ...
          'Interpreter','latex');
    
    colormap(ax, parula);
    cb = colorbar(ax);
    cb.Label.String = '$\Delta$';
    cb.Label.Interpreter = 'latex';
    
    view(135,25);
    box on;
    
    % opzionale: forza la scala colore per evidenziare che Δ < 0
    caxis([min(dvals) 0]);
    
    % salvataggio figura per LaTeX
    set(fig,'PaperPositionMode','auto');
    print(fig,'Delta_cubico_3D.png','-dpng','-r400');
