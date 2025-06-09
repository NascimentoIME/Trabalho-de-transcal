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

% ------------------ DOMÍNIO TEMPORAL ------------------
Fo_vals = logspace(-6, 2, 300);        % Número de Fourier (adimensional)
t_vals = Fo_vals * (L^2 / alpha);      % Tempo dimensional [s]

% ------------------ AUTOVALORES E COEFICIENTES (SERIE DE FOURIER) ------------------
N_values = [1, 5, 10, 50, 100];        % Número de termos
maxN = max(N_values);

% Equação transcendental: mu * cot(mu) + Bi = 0
f = @(mu) mu .* cot(mu) + Bi;

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

% ------------------ CÁLCULO DA TEMPERATURA MÉDIA ------------------
TempMedia_adim_all = zeros(length(Fo_vals), length(N_values));
TempMedia_dim_all = zeros(length(Fo_vals), length(N_values));

for p = 1:length(N_values)
    N = N_values(p);
    theta_med = zeros(size(Fo_vals));

    for k = 1:N
        theta_med = theta_med + ...
            2 * (mu(k)^2 + Bi^2) * (1 - cos(mu(k)))^2 .* exp(-mu(k)^2 * Fo_vals) ...
            / ((mu(k)^2 + Bi^2 + Bi) * mu(k)^2);
    end

    TempMedia_adim_all(:, p) = theta_med;
    TempMedia_dim_all(:, p) = Tinf + theta_med * (Tp - Tinf);
end

% ------------------ PLOTAGEM TEMPERATURA MÉDIA DIMENSIONAL ------------------
figure;
h = semilogx(t_vals, TempMedia_dim_all, 'LineWidth', 2); hold on;

% Curvas de +10% e -10% da curva com N = 100
T_ref = TempMedia_dim_all(:, end);          % curva com N = 100
T_plus10 = Tinf + 1.10 * (T_ref - Tinf);
T_minus10 = Tinf + 0.90 * (T_ref - Tinf);

h_plus = semilogx(t_vals, T_plus10, 'r--', 'LineWidth', 1.5);
h_minus = semilogx(t_vals, T_minus10, 'r--', 'LineWidth', 1.5);

% Legendas
leg_labels = arrayfun(@(n) sprintf('N = %d', n), N_values, 'UniformOutput', false);
leg_labels{end+1} = '+10%% de N = 100';
leg_labels{end+1} = '-10%% de N = 100';
legend([h; h_plus; h_minus], leg_labels, 'Location', 'Best', 'FontSize', 10);

xlabel('Tempo [s]', 'FontSize', 12);
ylabel('Temperatura média \bar{T}(t) [K]', 'FontSize', 12);
title('Convergência da Temperatura Média (Dimensional) com Diferentes N', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);