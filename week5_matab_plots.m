% Common Parameters
mu = 0.0035; % Viscosity of blood (Pa.s or g/(cm.s))
rho = 1.05; % Density of blood (g/cm^3)
L = 10; % Length of the vessel (cm)
dz = 1; % Differential element in z-direction (cm)

% Time settings
t = linspace(0, 10, 100); % Time (s)
dt = t(2) - t(1); % Time step

% Parameters for Blood Flow Rate
S_values = [0.10, 0.15, 0.20, 0.25, 0.30]; % Cross-sectional areas (cm^2)
Q = zeros(length(S_values), length(t)); % Initialize Q matrix
dPdz = 5; % Pressure gradient (mmHg/cm)

% Initial conditions
Q(:, 1) = 0; % Assuming initial flow rate is zero

% Solve differential equation for each S
figure;
hold on;
for j = 1:length(S_values)
    S = S_values(j);
    for i = 2:length(t)
        dQdt = -(4 * pi * mu / S) * Q(j, i-1) - (S / (2 * rho)) * dPdz;
        Q(j, i) = Q(j, i-1) + dQdt * dt; % Euler method
    end
    plot(t, Q(j, :), 'DisplayName', sprintf('S = %.2f cm^2', S));
end
hold off;
title('Blood Flow Rate versus Time for Different Cross-Sectional Areas');
xlabel('Time (s)');
ylabel('Blood Flow Rate (Q)');
legend;
grid on;

% Vessel radius values
R = linspace(0.1, 1, 100); % Radius of the vessel (cm)
Q_poiseuille = zeros(size(R)); % Initialize Q array

% Solve using Poiseuille’s equation
for i = 1:length(R)
    Q_poiseuille(i) = (pi * (R(i)^4) * dPdz) / (8 * mu * L); % Poiseuille's law
end

% Plot
figure;
plot(R, Q_poiseuille);
title('Figure 2: Variation of Blood Flow Rate with Vessel Radius Using Poiseuille’s Equation');
xlabel('Vessel Radius (cm)');
ylabel('Blood Flow Rate (Q)');
grid on;

% Different pressure gradients
pressure_gradients = [5, 10, 15, 20]; % Pressure gradient (mmHg/cm)
Q_gradients = zeros(length(pressure_gradients), length(t));

% Solve differential equation for each pressure gradient
figure;
hold on;
for j = 1:length(pressure_gradients)
    dPdz = pressure_gradients(j);
    Q_gradients(j, 1) = 0; % Initial condition
    for i = 2:length(t)
        dQdt = -(4 * pi * mu / S_values(3)) * Q_gradients(j, i-1) - (S_values(3) / (2 * rho)) * dPdz; % Example with S = 0.2 cm^2
        Q_gradients(j, i) = Q_gradients(j, i-1) + dQdt * dt;
    end
    plot(t, Q_gradients(j, :), 'DisplayName', sprintf('dP/dz = %d mmHg/cm', dPdz));
end
hold off;
title('Figure 3: Blood Flow Rate for Different Pressure Gradients');
xlabel('Time (s)');
ylabel('Blood Flow Rate (Q)');
legend;
grid on;


% Parameters for Blood Pressure
S_values = [0.10, 0.15, 0.20, 0.25, 0.30]; % Cross-sectional areas (cm^2)
P = zeros(length(S_values), length(t)); % Initialize P matrix

% Initial conditions
P(:, 1) = 100; % Assuming initial pressure is 100 mmHg

% Solve differential equation for each S
figure;
hold on;
for j = 1:length(S_values)
    S = S_values(j);
    for i = 2:length(t)
        dPdt = -(4 * mu / S^2) * P(j, i-1) - (4 * L * mu / (rho * S^2)) * dPdz;
        P(j, i) = P(j, i-1) + dPdt * dt; % Euler method
    end
    plot(t, P(j, :), 'DisplayName', sprintf('S = %.2f cm^2', S));
end
hold off;
title('Figure 4: Blood Pressure for Different Cross-Sectional Areas');
xlabel('Time (s)');
ylabel('Blood Pressure (P)');
legend;
grid on;

% Parameters
Q = 5; % Blood flow rate (cm^3/s)
mu = 0.0035; % Viscosity of blood (Pa.s or g/(cm.s))
L = 10; % Length of the vessel (cm)
R = linspace(0.1, 1, 100); % Radius of the vessel (cm)

% Poiseuille's equation for pressure drop
P = (8 * Q * mu * L) ./ (pi * R.^4);

% Plotting
figure;
plot(R, P, 'LineWidth', 2);
title('Figure 5: Variation of Blood Pressure with Vessel Radius using Poiseuille’s ');
xlabel('Vessel Radius (cm)');
ylabel('Blood Pressure Drop (ΔP)');
grid on;


% Parameters for Figure 7
Q = 5; % Blood flow rate (cm^3/s)
mu = 0.0035; % Viscosity of blood (Pa.s or g/(cm.s))
R = 0.2; % Fixed radius of the vessel (cm)
L = linspace(0, 10, 100); % Vessel length (cm), starting from 0 to 10

% Poiseuille's equation for pressure drop
P = (8 * Q * mu * L) / (pi * R^4);

% Slope calculation for proportionality
slope = (8 * mu) / (pi * R^4);

% Plotting
figure;
plot(L, P, 'LineWidth', 2);
title('Variation of Blood Pressure with Vessel Length using Poiseuille’s Equation');
xlabel('Vessel Length (cm)');
ylabel('Blood Pressure Drop (ΔP)');
grid on;
