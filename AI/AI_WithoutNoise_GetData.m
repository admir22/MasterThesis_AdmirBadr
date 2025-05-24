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

% Initializing the table to store results
num_cases = length(Te_values);
results = zeros(num_cases, 7); % 7 columns: Te, Nf, Speed, ia, Iq_2nd, Iq_4th, Rf

% Loop over all cases
for k = 1:num_cases
    % Load the current conditions
    Te = Te_values(k);
    Nf = Nf_values(k);
    speed_i = Speed_values(k);
    ia = ia_values(k);  % Added the new column `ia`
    Rf = Rf_values_OL(k); 

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

    out = sim('faulty_motor_model');

    time = tf/2:1/fs:tf-1/fs;
    i_d = interp1(out.simout1.time, out.simout1.signals.values(:,1),time, "linear");
    i_q = interp1(out.simout1.time, out.simout1.signals.values(:,2),time, "linear");

    % FFT parameters
    T = time(end) - time(1);
    N = length(time);
    Fs = (N-1)/T;
    fh = (-N/2:N/2-1)*(Fs/N); % Frequency vector

    % Fourier Transform
    I_d_fft = fft(i_d);  
    I_q_fft = fft(i_q);  

    % Extracting positive frequencies
    m = length(fh);
    f = fh(m/2+1:end);
    f = f(2:end); % Removing DC component

    % Positive spectrum
    m = length(I_q_fft);
    positiveiqfft = I_q_fft(2:m/2);

    % Identifying the 2nd harmonic
    [~, harmonique_idx] = max(abs(positiveiqfft)); 

    % Index for the 4th harmonic
    harmonique_idx4 = 2 * harmonique_idx;  

    % Checking indices to prevent overflow
    I_q_2nd = I_q_fft(1 + harmonique_idx);
    I_q_4th = I_q_fft(1 + harmonique_idx4);

    % Calculating normalized magnitudes
    Iq_2nd = 20 * log10(abs(I_q_2nd) / N);
    Iq_4th = 20 * log10(abs(I_q_4th) / N);

    % Storing the results
    results(k, :) = [Te, Nf, speed_i, ia, Iq_2nd, Iq_4th, Rf];

    % Displaying results for each case
    fprintf('Case %d: Te = %.2f, Nf = %d, Speed = %d, ia = %.2f, Rf = %.6f, Iq_2nd = %.5f, Iq_4th = %.5f\n', ...
        k, Te, Nf, speed_i, ia, Rf, Iq_2nd, Iq_4th);
end

% Final display as a table
disp('Results table:');
disp(array2table(results, ...
    'VariableNames', {'Te', 'Nf', 'Speed', 'ia', 'Iq_2nd', 'Iq_4th', 'Rf'}));

% Saving results for AI training
save('harmonics_data.mat', 'results');
