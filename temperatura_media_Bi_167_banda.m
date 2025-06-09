clc; clear; close all;

% Parâmetros fixos
Bi = 1.67;
Fo_vals = logspace(-6, 2 , 300); % ESCALA LOG no tempo
N_values = [1, 5, 10, 50, 100];  % Número de termos

% Função transcendental para autovalores
f = @(mu) mu .* cot(mu) + Bi;

% Pré-calcular autovalores (mu) e coeficientes (Ai)
maxN = max(N_values);
mu = zeros(1, maxN);
Ai = zeros(1, maxN);

for i = 1:maxN
    a = (i-1)*pi + 0.001;
    b = i*pi - 0.001;
    mu(i) = fzero(f, [a, b]);
end

for i = 1:maxN
    num = 2 * (mu(i)^2 + Bi^2) * (1 - cos(mu(i)));
    den = (mu(i)^2 + Bi^2 + Bi) * mu(i);
    Ai(i) = num / den;
end

% Matriz para armazenar a temperatura média para cada N
TempMedia_all = zeros(length(Fo_vals), length(N_values));

% Calcular temperatura média para cada N
for p = 1:length(N_values)
    N = N_values(p);
    TempMedia = zeros(size(Fo_vals));

    for k = 1:N
        TempMedia = TempMedia + ...
            2 * (mu(k)^2 + Bi^2) * (1 - cos(mu(k)))^2 .* exp(-mu(k)^2 * Fo_vals) ...
            / ((mu(k)^2 + Bi^2 + Bi) * mu(k)^2);
    end

    TempMedia_all(:, p) = TempMedia;
end

% Plotar todas as curvas em um único gráfico
figure;
h = semilogx(Fo_vals, TempMedia_all, 'LineWidth', 2); hold on;

% Curvas de +10% e -10% da curva com N = 100
theta_ref = TempMedia_all(:, end);          % curva com N = 100
theta_plus10 = 1.10 * theta_ref;
theta_minus10 = 0.90 * theta_ref;

h_plus = semilogx(Fo_vals, theta_plus10, 'r--', 'LineWidth', 1.5);
h_minus = semilogx(Fo_vals, theta_minus10, 'r--', 'LineWidth', 1.5);

% Legenda com nomes apropriados
leg_labels = arrayfun(@(n) sprintf('N = %d', n), N_values, 'UniformOutput', false);
leg_labels{end+1} = '+10% de N = 100';
leg_labels{end+1} = '-10% de N = 100';
legend([h; h_plus; h_minus], leg_labels, 'Location', 'Best', 'FontSize', 10);

xlabel('Fo (Número de Fourier)', 'FontSize', 12);
ylabel('\theta_{média}(Fo)', 'FontSize', 12);
title('Convergência da Temperatura Média com Diferentes N', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
