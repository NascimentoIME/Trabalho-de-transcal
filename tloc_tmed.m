clc; clear; close all;

% --- Parâmetros fixos ---
Bi = 1.67;
lambda = 0.6;                        % Posição fixa para temperatura local
Fo_vals = logspace(-6, 2, 300);      % Escala log no tempo
N = 100;                             % Número de termos

% --- Função transcendental para autovalores ---
f = @(mu) mu .* cot(mu) + Bi;

% --- Encontrar os autovalores mu ---
mu = zeros(1, N);
for i = 1:N
    a = (i - 1) * pi + 0.001;
    b = i * pi - 0.001;
    mu(i) = fzero(f, [a, b]);
end

% --- Calcular os coeficientes Ai para temperatura média ---
Ai = zeros(1, N);
for i = 1:N
    num = 2 * (mu(i)^2 + Bi^2) * (1 - cos(mu(i)));
    den = (mu(i)^2 + Bi^2 + Bi) * mu(i);
    Ai(i) = num / den;
end

% --- Calcular temperatura média θ_m(Fo) ---
TempMedia = zeros(1, length(Fo_vals));
for i = 1:N
    TempMedia = TempMedia + ...
        2 * (mu(i)^2 + Bi^2) * (1 - cos(mu(i)))^2 .* exp(-mu(i)^2 * Fo_vals) ...
        / ((mu(i)^2 + Bi^2 + Bi) * mu(i)^2);
end

% --- Calcular temperatura local θ(λ, Fo) ---
Ai_local = Ai;  % Pode reaproveitar Ai
TempLocal = zeros(1, length(Fo_vals));
for i = 1:N
    TempLocal = TempLocal + ...
        Ai_local(i) * sin(mu(i) * lambda) .* exp(-mu(i)^2 * Fo_vals);
end

% --- Plotar as curvas ---
figure;
semilogx(Fo_vals, TempLocal, 'k-', 'LineWidth', 2); hold on;
semilogx(Fo_vals, TempMedia, 'g-', 'LineWidth', 2);
xlabel('Fo (Número de Fourier)', 'FontSize', 12);
ylabel('\theta (Temperatura adimensional)', 'FontSize', 12);
title(sprintf('Comparação: Temperatura Local (λ = %.2f) vs. Temperatura Média', lambda), 'FontSize', 14);
legend({'Temperatura Local', 'Temperatura Média'}, 'Location', 'best');
grid on;
set(gca, 'FontSize', 12);
