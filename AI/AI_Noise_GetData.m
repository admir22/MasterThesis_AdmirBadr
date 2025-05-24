clear; clc; close all;

% Data for different cases
Te_values = [0, 0, 9.65, 9.65];
Rf = 0;
Speed_values = [500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000];
ia_values = [0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5];

Rf_values_OL = [0.095298357, 0.015888782, 0.001034171, 0.201119266, 0.04235185, 0.010600066, ...
             0.476688449, 0.079492148, 0.001034171, 1.005790838, 0.211420157, 0.051452884, ...
             0.953762807, 0.157572236, 0.001034171, 2.011472209, 0.419077365, 0.091013313, ...
             0.149183488, 0.029302598, 0.005326225, 0.290095747, 0.064450959, 0.019367738, ...
             0.7447346, 0.145043992, 0.024522019, 1.447639222, 0.318975078, 0.091970699, ...
             1.486039207, 0.284669681, 0.039594496, 2.887719507, 0.625906837, 0.163029848];

% Initializing the data for storing results
results = [];

% Loop over each case
for i = 1:length(Te_values)
    Te = Te_values(i);
    Nf = Nf_values(i);
    speed_i = Speed_values(i);
    ia = ia_values(i);
    Rf = Rf_values_OL(i); 

    % Initializing parameters
    mu = (Nf/378); % Fault severity
    
    % Simulation
    faulty_PMSM_load;
    vdq;
    
    seedA = 0;
    seedB = 1;
    seedC = 2;
    seedW = 3;
    seedId = 4;
    seedIq = 5;
    seedCL = 6;
    
    out = sim('faulty_motor_model');
    
    % Simulation parameters
    time = 0.2:1/fs:tf-1/fs;
    T_fundamental = 60 / (speed_i * 3);
    N_periods = floor((time(end)-time(1)) / T_fundamental);
    
    % Calculating instantaneous harmonic values
    for k = 1:N_periods
        t_start = (k - 1) * T_fundamental + 0.2;
        t_end = k * T_fundamental + 0.2;
        interval = t_start:1/fs:t_end - 1/fs;
        
        i_q = interp1(out.noiseharmo.time, out.noiseharmo.signals.values(:,2), interval, "linear");
        
        % Fourier Transform
        I_q_fft = fft(i_q);
        
        % Identifying the 2nd and 4th harmonics
        positiveiqfft = I_q_fft(2:end/2);
        [~, harmonique_idx] = max(abs(positiveiqfft));
        harmonique_idx4 = 2 * harmonique_idx;
        
        Iq_2nd_inst = 20 * log10(abs(I_q_fft(1 + harmonique_idx)) / length(interval));
        Iq_4th_inst = 20 * log10(abs(I_q_fft(1 + harmonique_idx4)) / length(interval));
        
        % Storing each instantaneous value as a separate point
        results = [results; Te, Nf, speed_i, ia, Iq_2nd_inst, Iq_4th_inst];
    end
end

% Convert to table for display and saving
dataset = array2table(results, 'VariableNames', {'Te', 'Nf', 'Speed', 'ia', 'Iq_2nd_inst', 'Iq_4th_inst'});

AD = table2array(dataset);

% Saving the results
save('harmonics_data_instantaneousAD.mat', 'AD');

disp('Instantaneous harmonic data successfully saved.');

% Displaying the first few rows for verification
disp(dataset);
