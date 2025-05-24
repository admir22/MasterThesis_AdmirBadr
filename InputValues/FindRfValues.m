clear; clc; close all;

% Parameter definition
Te_values = [0, 9.65]; % Electromagnetic torque values (Nm)
mu_values = [1/378, 5/378, 10/378]; % Values of the mutual inductance to phase inductance ratio
speed_values = [500, 1000]; % Speed values (rpm)
If_targets = [0.5, 2, 5]; % Target flux values (A)

% Storing the results
results = [];

for Te = Te_values
    for mu = mu_values
        for speed_i = speed_values
            fprintf('Simulating for Te=%.2f, mu=%.6f, speed=%d\n', Te, mu, speed_i);
            
            % Load faulty PMSM model and calculate vdq
            faulty_PMSM_load;
            vdq;
            
            Rf_values = zeros(size(If_targets));
            
            for i = 1:length(If_targets)
                target_If = If_targets(i);
                
                % Finding Rf using fminbnd to bound the interval
                try
                    Rf_values(i) = fminbnd(@(Rf) abs(simulate_motor(Rf, Te, mu, speed_i) - target_If), 0.001, 7);
                catch
                    Rf_values(i) = NaN; % In case of an error, assign NaN
                    warning('Convergence issue for Te=%.2f, mu=%.6f, speed=%d, If_target=%.2f', Te, mu, speed_i, target_If);
                end
            end
            
            results = [results; repmat([Te, mu, speed_i], length(If_targets), 1), If_targets', Rf_values'];
        end
    end
end

disp('Simulation results:');
disp(array2table(results, 'VariableNames', {'Te', 'mu', 'speed', 'If_target', 'Rf'}));

% Save results to an Excel file
filename = 'Rf_results.xlsx';
writematrix(results, filename);
fprintf('Results saved in %s\n', filename);

% Motor simulation function
function If_value = simulate_motor(Rf, Te, mu, speed_i)
    assignin('base', 'Rf', Rf);
    assignin('base', 'Te', Te);
    assignin('base', 'mu', mu);
    assignin('base', 'speed_i', speed_i);
    
    try
        out = sim('faulty_motor_model');
        If_value = out.simout2.signals.values(end);
    catch
        If_value = NaN; % If the simulation fails, return NaN
    end
end
