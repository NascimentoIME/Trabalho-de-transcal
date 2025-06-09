
%%%% tijolo_comum_06_semifinal_transcal


clc; clear; close all;

% ------------------ PROPRIEDADES DO MATERIAL (TIJOLO COMUM) ------------------
k = 0.72;               % condutividade térmica [W/m·K]
rho = 1920;             % densidade [kg/m³]
cp = 835;               % calor específico [J/kg·K]
L = 0.09;               % espessura da parede [m]
alpha = k / (rho * cp); % difusividade térmica [m²/s]

% ------------------ CONDIÇÕES DE CONTORNO - CONVECÇÃO NATURAL ------------------
h = 3.0357;             % coeficiente de convecção [W/m²·K]
Bi = h * L / k;         % número de Biot (adimensional)

% ------------------ CONDIÇÕES INICIAIS E DE REFERÊNCIA ------------------
Tp = 300;               % Temperatura inicial do tijolo [K]
Tinf = 290;             % Temperatura do ar ambiente [K]

% ------------------ DOMÍNIO TEMPORAL E POSICIONAL ------------------
Fo_values = logspace(-6, 2, 100);    % Número de Fourier (escala logarítmica)
t_values = Fo_values * (L^2 / alpha); % Tempo dimensional [s]

lambda = 0.6;                      % Posição adimensional x/L
x_real = lambda * L;              % Posição física real [m]

% ------------------ TERMOS DA SÉRIE DE FOURIER ------------------
N_values = [1, 10, 20, 40, 100];     % Número de termos para análise de convergência
T_vals = zeros(length(Fo_values), length(N_values)); % Temperatura dimensional [K]

for k_idx = 1:length(Fo_values)
    Fo = Fo_values(k_idx);

    for j = 1:length(N_values)
        N = N_values(j);

        % Resolver equação transcendental para obter autovalores mu
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

        % Solução adimensional theta (lambda = x/L fixado)
        theta = sum(A .* sin(mu * lambda) .* exp(-mu.^2 * Fo));

        % Solução dimensional
        T_vals(k_idx, j) = theta * (Tp - Tinf) + Tinf;
    end
end

% ------------------ GRÁFICO 1: Temperatura dimensional T(x,t) ------------------
figure;
semilogx(t_values, T_vals, 'LineWidth', 2);
xlabel('t (s)', 'FontSize', 12);
ylabel(sprintf('Temperatura T(x=%.3f m, t) [K]', x_real), 'FontSize', 12);
title(sprintf('Convergência da Solução em x = %.3f m com Diferentes N', x_real), 'FontSize', 14);

legend(arrayfun(@(n) sprintf('N = %d', n), N_values, 'UniformOutput', false), ...
       'Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);

% ------------------ GRÁFICO 2: Erro relativo em relação a N=100 ------------------
T_ref = T_vals(:, end);         % Referência (N=100)
T_ref(T_ref == 0) = 1e-12;      % Evita divisão por zero
erro_rel = abs(T_vals - T_ref) ./ abs(T_ref); % Erro relativo

figure;
semilogx(t_values, erro_rel(:, 1:end-1)*100, 'LineWidth', 2);
xlabel('t (s)', 'FontSize', 12);
ylabel('Erro relativo (%)', 'FontSize', 12);
title(sprintf('Erro Relativo em x = %.3f m em Relação a N = 100 Termos', x_real), 'FontSize', 14);

legend(arrayfun(@(n) sprintf('N = %d', n), N_values(1:end-1), 'UniformOutput', false), ...
       'Location', 'Best', 'FontSize', 10);
ylim([0 10]);
grid on;
set(gca, 'FontSize', 12);

% ------------------ GRÁFICO 3: Faixa de tolerância ±10% ------------------
T_plus10 = Tinf + 1.10 * (T_ref - Tinf);
T_minus10 = Tinf + 0.90 * (T_ref - Tinf);

figure;
semilogx(t_values, T_ref, 'k-', 'LineWidth', 2); hold on;
semilogx(t_values, T_plus10, 'r--', 'LineWidth', 1.5);
semilogx(t_values, T_minus10, 'r--', 'LineWidth', 1.5);
colors = lines(length(N_values)-1);

for j = 1:length(N_values)-1
    semilogx(t_values, T_vals(:, j), '-', 'Color', colors(j,:), 'LineWidth', 1.2);
end

xlabel('t (s)', 'FontSize', 12);
ylabel(sprintf('Temperatura T(x=%.3f m, t) [K]', x_real), 'FontSize', 12);
title(sprintf('Curvas de T(x=%.3f m, t) e Faixas de Tolerância (±10%%)', x_real), 'FontSize', 14);

legend({'N=100 (referência)', '+10%', '-10%', 'N=1', 'N=10', 'N=20', 'N=40'}, ...
    'Location', 'best');
grid on;
set(gca, 'FontSize', 12);
