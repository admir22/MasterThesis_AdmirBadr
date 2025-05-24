clc;

% System data
R = Rs; % Phase resistance (Ohms)
Ld = Ld_0; % Direct inductance (Henries)
Lq = Lq_0; % Quadrature inductance (Henries)
Psi_PM = psi_pm; % Permanent flux (Webers)
p = 3; % Number of pole pairs
Te_target = 9.65; % Desired electromagnetic torque (Nm)
w = speed_i * 2 * pi / 60; % Angular speed (rad/s)
i_max = 2.7*sqrt(2); % Maximum allowable current (A)

% MTPA strategy calculation
% Total reference current
i_m_ref = (2 * Te_target) / (3 * p * Psi_PM);

% Limiting the reference current
i_m = min(i_m_ref, i_max);

% Calculating i_d according to MTPA
i_d_mtpa = (Psi_PM / (4 * (Lq - Ld))) - sqrt((Psi_PM / (4 * (Lq - Ld)))^2 + (i_m^2 / 2));

% Calculating i_q according to MTPA
i_q_mtpa = sqrt(i_m^2 - i_d_mtpa^2);

% Calculating vd and vq voltages using i_d_mtpa and i_q_mtpa
syms vd vq real

% Static equations
eq1 = vd == R * i_d_mtpa - p * w * Lq * i_q_mtpa;
eq2 = vq == R * i_q_mtpa + p * w * Ld * i_d_mtpa + p * w * Psi_PM;

% Solving for vd and vq
[sol_vd, sol_vq] = solve([eq1, eq2], [vd, vq]);

% Converting to numeric values
vd_value = double(vpa(sol_vd, 6));
vq_value = double(vpa(sol_vq, 6));

% Display results as integers or usable numbers
disp('Solutions ready for Simulink:');
disp(['vd (MTPA) = ', num2str(vd_value)]);
disp(['vq (MTPA) = ', num2str(vq_value)]);
disp(['id (MTPA) = ', num2str(i_d_mtpa)]);
disp(['iq (MTPA) = ', num2str(i_q_mtpa)]);
