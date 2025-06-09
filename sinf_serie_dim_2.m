clc; clear; close all;

% ------------------ PROPRIEDADES DO MATERIAL (TIJOLO COMUM) ------------------
k = 0.72; rho = 1920; cp = 835;
L = 0.09;                                % espessura da parede [m]
alpha = k / (rho * cp);                 % difusividade térmica [m²/s]
h = 3.0357;                              % coef. de convecção [W/m²·K]
Bi = h * L / k;                          % número de Biot

% ------------------ CONDIÇÕES INICIAIS E CONTORNO ------------------
Tp = 300; Tinf = 290;

% ------------------ DOMÍNIO TEMPORAL ------------------
Fo_vals = logspace(-6, 2, 300);
t_vals = Fo_vals * (L^2 / alpha);

% ------------------ POSIÇÃO FIXA (x = L) ------------------
lambda = 1.0; 
x_real = lambda * L;  % = 0.09 m

% ------------------ 1. SOLUÇÃO ANALÍTICA (SEMI-INFINITA) COM CONVECÇÃO ------------------
x_espelhado = L - x_real;  % espelhamento para x = L

term1 = erfc(x_espelhado ./ (2 * sqrt(alpha * t_vals)));
exp_term = exp((h * x_espelhado)/k + ((h^2) * alpha .* t_vals) / k^2);
term2 = erfc(x_espelhado ./ (2 * sqrt(alpha * t_vals)) + (h * sqrt(alpha * t_vals)) / k);

theta_conv = term1 - exp_term .* term2;
T_convectiva = Tp + theta_conv * (Tinf - Tp);  % temperatura com convecção na face oposta

% ------------------ 2. SOLUÇÃO DA SÉRIE DE FOURIER (N = 100) ------------------
N = 100;
T_serie = zeros(size(t_vals));

for k_idx = 1:length(Fo_vals)
    Fo = Fo_vals(k_idx);
    
    % Encontrar autovalores mu (raízes da equação transcendental)
    f = @(mu) mu .* cot(mu) + Bi;
    mu = zeros(1, N);
    for i = 1:N
        a = (i - 1) * pi + 0.001;
        b = i * pi - 0.001;
        mu(i) = fzero(f, [a, b]);
    end

    % Coeficientes A_n
    A = zeros(1, N);
    for i = 1:N
        num = 2 * (mu(i)^2 + Bi^2) * (1 - cos(mu(i)));
        den = (mu(i)^2 + Bi^2 + Bi) * mu(i);
        A(i) = num / den;
    end

    % Solução adimensional e dimensional
    theta = sum(A .* sin(mu * lambda) .* exp(-mu.^2 * Fo));
    T_serie(k_idx) = theta * (Tp - Tinf) + Tinf;
end

% ------------------ PLOTAGEM COMPARATIVA ------------------
figure;
semilogx(t_vals, T_serie, 'r-', 'LineWidth', 2); hold on;
semilogx(t_vals, T_convectiva, 'b--', 'LineWidth', 2);
grid on;

xlabel('Tempo [s]', 'FontSize', 12);
ylabel(sprintf('Temperatura T(x=%.3f m, t) [K]', x_real), 'FontSize', 12);
title(sprintf('Comparação entre Série de Fourier (N=100) e Solução Semi-infinita (espelhada) em x = %.3f m', x_real), 'FontSize', 14);
legend('Série de Fourier (N=100)', 'Semi-infinita com Convecção', 'Location', 'Best');
set(gca, 'FontSize', 12);
