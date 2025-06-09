clc; clear; close all;

% ------------------ PROPRIEDADES DO MATERIAL (TIJOLO COMUM) ------------------
k = 0.72;               % condutividade térmica [W/m·K]
rho = 1920;             % densidade [kg/m³]
cp = 835;               % calor específico [J/kg·K]
L = 0.09;               % espessura da parede [m]
alpha = k / (rho * cp); % difusividade térmica [m²/s]

% ------------------ CONDIÇÕES INICIAIS E DE CONTORNO ------------------
Tp = 300;               % Temperatura inicial do tijolo [K]
Tinf = 290;             % Temperatura do meio externo [K]

% ------------------ CONDIÇÕES DE CONTORNO - CONVECÇÃO NATURAL ------------------
h = 3.0357;             % coeficiente de convecção [W/m²·K]
Bi = h * L / k;         % número de Biot (adimensional)

% ------------------ POSIÇÃO FIXA ------------------
lambda = 0.1;           % posição adimensional x/L
x_real = lambda * L;    % posição física real [m]

% ------------------ DOMÍNIO TEMPORAL ------------------
Fo_vals = logspace(-6, 2, 300);       % Número de Fourier (adimensional)
t_vals = Fo_vals * (L^2 / alpha);     % Tempo dimensional [s]

% ------------------ SOLUÇÃO SEMI-INFINITA ------------------
theta_semi = erf( x_real ./ (2 * sqrt(alpha * t_vals)) );   % solução adimensional semi-infinita
T_semi = Tinf + theta_semi * (Tp - Tinf);                   % solução dimensional semi-infinita

% ------------------ SOLUÇÃO EM SÉRIE (N=100 termos) ------------------
N = 100;                    % número de termos da série
T_series = zeros(size(t_vals));

% Função transcendental para autovalores mu
f = @(mu) mu .* cot(mu) + Bi;

% Encontrar autovalores mu_i
mu = zeros(1, N);
for i = 1:N
    a = (i - 1)*pi + 0.001;
    b = i*pi - 0.001;
    mu(i) = fzero(f, [a b]);
end

% Calcular coeficientes A_n
A = zeros(1, N);
for i = 1:N
    num = 2 * (mu(i)^2 + Bi^2) * (1 - cos(mu(i)));
    den = (mu(i)^2 + Bi^2 + Bi) * mu(i);
    A(i) = num / den;
end

% Calcular temperatura para cada tempo t_vals
for k_idx = 1:length(t_vals)
    Fo = Fo_vals(k_idx);
    theta_series = sum( A .* sin(mu * lambda) .* exp(-mu.^2 * Fo) );
    T_series(k_idx) = Tinf + theta_series * (Tp - Tinf);
end

% ------------------ PLOTAGEM ------------------
figure;
semilogx(t_vals, T_semi, 'r-', 'LineWidth', 2); hold on;
semilogx(t_vals, T_series, 'b-', 'LineWidth', 2);
grid on;

xlabel('Tempo [s]', 'FontSize', 12);
ylabel(sprintf('Temperatura T(x=%.3f m, t) [K]', x_real), 'FontSize', 12);
title('Comparação: Solução Semi-infinita vs Série com 100 termos', 'FontSize', 14);
legend({'Solução Semi-infinita (erf)', 'Solução em Série (N=100)'}, 'Location', 'Best');
set(gca, 'FontSize', 12);
