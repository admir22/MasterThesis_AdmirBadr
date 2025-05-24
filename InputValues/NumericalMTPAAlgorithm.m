% Define physical parameters

% Define the numerical values of known variables
R = Rs;                 % Stator resistance in Ohms
Psi_PM = psi_pm;        % Permanent magnet flux in Weber
Lq = Lq_0;              % Quadrature inductance in Henry
Ld = Ld_0;              % Direct inductance in Henry
p = 3;                  % Number of pole pairs
omega = speed_i * 2 * pi / 60; % Electrical angular speed in rad/s
Id_nom = -1.6;          % Nominal value of id in A
Iq_nom = 3.5;           % Nominal value of iq in A

% Objective function to minimize im²
objective = @(x) x(1)^2 + x(2)^2; % x(1) = id, x(2) = iq

% Initial guess for [id, iq, vq]
x0 = [0; 0; 0];

% Bounds for id and iq
lb = [Id_nom, -Iq_nom, -Inf]; % Lower bounds
ub = [-Id_nom, Iq_nom, Inf];  % Upper bounds

% Options for fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Numerical solution with constraints
[x_sol, fval, exitflag] = fmincon(@(x) objective(x), x0, [], [], [], [], lb, ub, @(x) constraints(x, Te, p, Psi_PM, Ld, Lq, R, omega), options);

% Check for convergence
if exitflag <= 0
    disp('The solution did not converge. Try a different initial guess or check the equations.');
else
    % Extract solutions
    id_sol = x_sol(1);
    iq_sol = x_sol(2);
    vq_sol = x_sol(3);

    % Calculate vd with the found values
    vd_sol = R * id_sol - p * omega * Lq * iq_sol;

    % Check the calculated torque
    Te_calc = (3/2) * p * (Psi_PM * iq_sol + (Ld - Lq) * id_sol * iq_sol);

    % Display results
    fprintf('Solution found:\n');
    fprintf('id = %.4f A\n', id_sol);
    fprintf('iq = %.4f A\n', iq_sol);
    fprintf('vd = %.4f V\n', vd_sol);
    fprintf('vq = %.4f V\n', vq_sol);
    fprintf('Calculated Te = %.4f Nm (target = %.4f Nm)\n', Te_calc, Te);
end

% Local function to define nonlinear constraints
function [c, ceq] = constraints(x, Te, p, Psi_PM, Ld, Lq, R, omega)
    id = x(1);
    iq = x(2);
    vq = x(3); % vq
    % Equality constraints
    ceq = [
        Te - (3/2) * p * (Psi_PM * iq + (Ld - Lq) * id * iq); % Torque equation
        vq - (R * iq + p * omega * Ld * id + p * omega * Psi_PM) % vq equation
    ];
    % No inequality constraints other than the bounds
    c = [];
end
