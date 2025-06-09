clc; clear; close all;

% ------------------ PROPRIEDADES DO MATERIAL (TIJOLO COMUM) ------------------
k = 0.72;               % [W/m·K]
rho = 1920;             % [kg/m³]
cp = 835;               % [J/kg·K]
L = 0.09;               % [m]
alpha = k / (rho * cp); % [m²/s]

% ------------------ CONDIÇÕES DE CONTORNO ------------------
h = 3.0357;             % [W/m²·K]
Bi = h * L / k;         % Número de Biot

% ------------------ CONDIÇÕES INICIAIS ------------------
Tp = 300;               % Temperatura inicial [K]
Tinf = 290;             % Temperatura ambiente [K]

% ------------------ DOMÍNIO TEMPORAL ------------------
Fo_vals = logspace(-6, 2, 300);        % Número de Fourier
t_vals = Fo_vals * (L^2 / alpha);      % Tempo real [s]

% ------------------ SOLUÇÃO: CAPACITÂNCIA GLOBAL ------------------
theta_cap = exp(-Bi * Fo_vals);        % Temperatura adimensional

% ------------------ SOLUÇÃO: SÉRIE COM N = 100 TERMOS ------------------
N = 100;
lambda = 0.9;                          % Posição adimensional
mu = zeros(1, N);                      % Autovalores
Ai = zeros(1, N);                      % Coeficientes

% Encontrar autovalores mu_n
f = @(mu) mu .* cot(mu) + Bi;
for i = 1:N
    a = (i - 1) * pi + 0.001;
    b = i * pi - 0.001;
    mu(i) = fzero(f, [a, b]);
end

% Coeficientes Ai
for i = 1:N
    num = 2 * (mu(i)^2 + Bi^2) * (1 - cos(mu(i)));
    den = (mu(i)^2 + Bi^2 + Bi) * mu(i);
    Ai(i) = num / den;
end

% Cálculo da temperatura adimensional via série
theta_series = zeros(size(Fo_vals));
for k = 1:length(Fo_vals)
    Fo = Fo_vals(k);
    theta_series(k) = sum(Ai .* sin(mu * lambda) .* exp(-mu.^2 * Fo));
end

% ------------------ PLOTAGEM COMPARATIVA ------------------
figure;
semilogx(t_vals, theta_cap, 'b--', 'LineWidth', 2); hold on;
semilogx(t_vals, theta_series, 'k-', 'LineWidth', 2);

xlabel('Tempo [s]', 'FontSize', 12);
ylabel('\theta = (T - T_{\infty}) / (T_p - T_{\infty})', 'FontSize', 14);
title('Capacitância Global vs Série (N = 100)', 'FontSize', 14);
legend('Capacitância Global', 'Série N = 100', 'Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);
