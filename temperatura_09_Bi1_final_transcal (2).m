clc; clear; close all;


%nome para salvar: temperatura_09_Bi1_final_transcal.m


% --- Parâmetros fixos ---
Bi = 5;
lambda = 0.9;
Fo_values = logspace(-6, 2, 100);    % Número de Fourier (escala log)
N_values = [1, 10, 20, 40, 100];    % Quantidade de termos da série

% --- Inicializar matriz para armazenar os valores de theta ---
theta_vals = zeros(length(Fo_values), length(N_values));

% --- Cálculo da solução em série para cada Fo e cada N ---
for k = 1:length(Fo_values)
    Fo = Fo_values(k);
    for j = 1:length(N_values)
        N = N_values(j);

        % Função transcendental para encontrar mu_n
        f = @(mu) mu .* cot(mu) + Bi;

        % Encontrar autovalores mu_n via fzero em cada intervalo (n*pi, (n+1)*pi)
        mu = zeros(1, N);
        for i = 1:N
            a = (i - 1) * pi + 0.001;
            b = i * pi - 0.001;
            mu(i) = fzero(f, [a, b]);
        end

        % Coeficientes A_n
        Ai = zeros(1, N);
        for i = 1:N
            num = 2 * (mu(i)^2 + Bi^2) * (1 - cos(mu(i)));
            den = (mu(i)^2 + Bi^2 + Bi) * mu(i);
            Ai(i) = num / den;
        end

        % Cálculo da temperatura adimensional theta
        theta = sum(Ai .* sin(mu * lambda) .* exp(-mu.^2 * Fo));
        theta_vals(k, j) = theta;
    end
end

% --- Cálculo do erro relativo em relação à referência N=100 ---
ref = theta_vals(:, end);
ref(ref == 0) = 1e-12; % Evitar divisão por zero
erro_rel = abs(theta_vals - ref) ./ abs(ref);

% --- GRÁFICO 1: Temperatura adimensional para diferentes N ---
figure;
semilogx(Fo_values, theta_vals, 'LineWidth', 2);
legend(arrayfun(@(n) sprintf('N = %d', n), N_values, 'UniformOutput', false), ...
       'Location', 'Best', 'FontSize', 10);
xlabel('Fo (Número de Fourier)', 'FontSize', 12);
ylabel('\theta(\lambda = 0.1, Fo)', 'FontSize', 12);
title('Convergência da Solução com Diferentes N', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);

% --- GRÁFICO 2: Erro relativo (%) em relação a N=100 ---
figure;
semilogx(Fo_values, erro_rel(:, 1:end-1)*100, 'LineWidth', 2);
legend(arrayfun(@(n) sprintf('N = %d', n), N_values(1:end-1), 'UniformOutput', false), ...
       'Location', 'Best', 'FontSize', 10);
xlabel('Fo (Número de Fourier)', 'FontSize', 12);
ylabel('Erro relativo (%)', 'FontSize', 12);
title('Erro Relativo em Relação à Solução com N = 100 Termos', 'FontSize', 14);
ylim([0 10]); % Limitar erro relativo para melhor visualização
grid on;
set(gca, 'FontSize', 12);

% Plotagem


theta_ref = theta_vals(:, end);  % N = 100

theta_plus10 = 1.10 * theta_ref;
theta_minus10 = 0.90 * theta_ref;

figure;
semilogx(Fo_values, theta_ref, 'k-', 'LineWidth', 2); hold on;
semilogx(Fo_values, theta_plus10, 'r--', 'LineWidth', 1.5);
semilogx(Fo_values, theta_minus10, 'r--', 'LineWidth', 1.5);
colors = lines(length(N_values)-1);

for j = 1:length(N_values)-1
    semilogx(Fo_values, theta_vals(:, j), '-', 'Color', colors(j,:), 'LineWidth', 1.2);
end

xlabel('Fo (Número de Fourier)', 'FontSize', 12);
ylabel('\theta (Temperatura adimensional)', 'FontSize', 12);
title('Curvas de \theta e Faixas de Tolerância (\pm10%)', 'FontSize', 14);
legend({'N=100 (referência)', '+10%', '-10%', 'N=1', 'N=10', 'N=20', 'N=40'}, ...
    'Location', 'best');
grid on;
set(gca, 'FontSize', 12);
