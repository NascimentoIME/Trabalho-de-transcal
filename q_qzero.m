clc; clear; close all;

% --- Parâmetros fixos ---
Bi = 1.67;
Fo_vals = logspace(-6, 2 , 300);    % Escala log no tempo
N_values = [1, 5, 10, 50, 100];     % Número de termos na série

% --- Função transcendental para autovalores ---
f = @(mu) mu .* cot(mu) + Bi;

% --- Pré-cálculo de autovalores (mu) e coeficientes (Ai) ---
maxN = max(N_values);
mu = zeros(1, maxN);
Ai = zeros(1, maxN);

for i = 1:maxN
    a = (i - 1) * pi + 0.001;
    b = i * pi - 0.001;
    mu(i) = fzero(f, [a, b]);
end

for i = 1:maxN
    num = 2 * (mu(i)^2 + Bi^2) * (1 - cos(mu(i)));
    den = (mu(i)^2 + Bi^2 + Bi) * mu(i);
    Ai(i) = num / den;
end

% --- Temperatura média e Q/Q0 ---
TempMedia_all = zeros(length(Fo_vals), length(N_values));
Q_Q0_all = zeros(length(Fo_vals), length(N_values));

for p = 1:length(N_values)
    N = N_values(p);
    TempMedia = zeros(size(Fo_vals));

    for k = 1:N
        TempMedia = TempMedia + ...
            2 * (mu(k)^2 + Bi^2) * (1 - cos(mu(k)))^2 .* exp(-mu(k)^2 * Fo_vals) ...
            / ((mu(k)^2 + Bi^2 + Bi) * mu(k)^2);
    end

    TempMedia_all(:, p) = TempMedia;
    Q_Q0_all(:, p) = 1 - TempMedia;  % Razão adimensional de calor
end

% --- Plotar Q/Q0 ---
figure;
semilogx(Fo_vals, Q_Q0_all, 'LineWidth', 2); hold on;

legend(arrayfun(@(n) sprintf('N = %d', n), N_values, 'UniformOutput', false), ...
       'Location', 'Best', 'FontSize', 10);

xlabel('Fo (Número de Fourier)', 'FontSize', 12);
ylabel('Q/Q_0', 'FontSize', 12);
title(' Q/Q_0 vs Fo', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
