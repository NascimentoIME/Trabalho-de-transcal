clc; clear; close all;

% ------------------ PROPRIEDADES DO MATERIAL (TIJOLO COMUM) ------------------
k = 0.72;               % condutividade térmica [W/m·K]
rho = 1920;             % densidade [kg/m³]
cp = 835;               % calor específico [J/kg·K]
L = 0.09;               % espessura da parede [m]
alpha = k / (rho * cp); % difusividade térmica [m²/s]

% ------------------ CONDIÇÕES DE CONTORNO - CONVECÇÃO NATURAL ------------------
h = 3.0357;             % coef. de convecção [W/m²·K]
Bi = h * L / k;         % número de Biot

% ------------------ TEMPERATURAS (para energia dimensional) ------------------
Tp = 300;               % Temperatura inicial [K]
Tinf = 290;             % Temperatura ambiente [K]
DeltaT = Tp - Tinf;     % Diferença de temperatura [K]

% Energia térmica inicial (por m² de parede)
Q0 = rho * cp * L * DeltaT;   % [J/m²]

% ------------------ Parâmetros da solução ------------------
Fo_vals = logspace(-6, 2 , 300);   % Número de Fourier (tempo adimensional)
tempo_segundos = Fo_vals * L^2 / alpha;  % Tempo real [s]
N_values = [1, 5, 10, 50, 100];    % Número de termos na série

% --- Função transcendental para autovalores ---
f = @(mu) mu .* cot(mu) + Bi;

% --- Pré-cálculo de autovalores e coeficientes ---
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

% --- Temperatura média adimensional e energia real Q(t) ---
TempMedia_all = zeros(length(Fo_vals), length(N_values));
Q_real_all = zeros(length(Fo_vals), length(N_values));   % Energia real (J/m²)

for p = 1:length(N_values)
    N = N_values(p);
    TempMedia = zeros(size(Fo_vals));

    for k = 1:N
        TempMedia = TempMedia + ...
            2 * (mu(k)^2 + Bi^2) * (1 - cos(mu(k)))^2 .* exp(-mu(k)^2 * Fo_vals) ...
            / ((mu(k)^2 + Bi^2 + Bi) * mu(k)^2);
    end

    TempMedia_all(:, p) = TempMedia;
    Q_real_all(:, p) = Q0 * (1 - TempMedia);  % Energia em J/m²
end

% --- Plotar energia térmica em função de tempo real ---
figure;
semilogx(tempo_segundos, Q_real_all, 'LineWidth', 2); hold on;

legend(arrayfun(@(n) sprintf('N = %d termos', n), N_values, 'UniformOutput', false), ...
       'Location', 'Best', 'FontSize', 10);

xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Energia Térmica Q(t) [J/m²]', 'FontSize', 12);
title('Energia Térmica Liberada pelo Sólido em Função do Tempo Real', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);

