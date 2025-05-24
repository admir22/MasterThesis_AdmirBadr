clear; clc; close all;

% Data from the table
Te_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
             9.65, 9.65, 9.65, 9.65, 9.65, 9.65, 9.65, 9.65, 9.65, 9.65, 9.65, 9.65, ...
             9.65, 9.65, 9.65, 9.65, 9.65, 9.65];
Nf_values = [1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, ...
             1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, ...
             10, 10, 10, 10, 10, 10];
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
Rf_values_CL = [0.095242124, 0.015864692, 0.001034171, 0.201071759, 0.042326903, 0.010565182, ...
             0.476038102, 0.07874583, 0.001034171, 1.004883506, 0.210301675, 0.05009258, ...
             0.951021885, 0.154219563, 0.001034171, 2.007495, 0.413986259, 0.083129123, ...
             0.149152878, 0.029265927, 0.00529274, 0.289947298, 0.064424072, 0.019294789, ...
             0.743692388, 0.143913688, 0.023190089, 1.445937047, 0.317072488, 0.089633134, ...
             1.481707035, 0.279600817, 0.032856456, 2.88060868, 0.617146849, 0.150810877];
Rf_values_OL = [0.095298357, 0.015888782, 0.001034171, 0.201119266, 0.04235185, 0.010600066, ...
             0.476688449, 0.079492148, 0.001034171, 1.005790838, 0.211420157, 0.051452884, ...
             0.953762807, 0.157572236, 0.001034171, 2.011472209, 0.419077365, 0.091013313, ...
             0.149183488, 0.029302598, 0.005326225, 0.290095747, 0.064450959, 0.019367738, ...
             0.7447346, 0.145043992, 0.024522019, 1.447639222, 0.318975078, 0.091970699, ...
             1.486039207, 0.284669681, 0.039594496, 2.887719507, 0.625906837, 0.163029848];
thresholdlist = [0.0160, 0.0160, 0.0160, 0.0213, 0.0213, 0.0213, ...
                 0.0160, 0.0160, 0.0160, 0.0213, 0.0213, 0.0213, ...
                 0.0160, 0.0160, 0.0160, 0.0213, 0.0213, 0.0213, ...
                 0.0160, 0.0160, 0.0160, 0.0213, 0.0213, 0.0213, ...
                 0.0160, 0.0160, 0.0160, 0.0213, 0.0213, 0.0213, ...
                 0.0160, 0.0160, 0.0160, 0.0213, 0.0213, 0.0213];

% Initializing results
results = [];

% Loop over each case given in the table
for i = 1:length(Te_values)
    Te = Te_values(i);
    Nf = Nf_values(i);
    speed_i = Speed_values(i);
    ia = ia_values(i);  % Added the new column `ia`
    Rf = Rf_values_OL(i); 
    threshold = thresholdlist(i);

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

    % Simulation parameters
    time = 0.2:1/fs:tf-1/fs; 
    T_fundamental = 60 / (speed_i * 3); % Fundamental period
    N_periods = floor((time(end)-time(1)) / T_fundamental);

    Iq = zeros(1, N_periods);
    Id = zeros(1, N_periods);

    % Calculating unbalance percentages
    for k = 1:N_periods
        t_start = (k - 1) * T_fundamental + 0.2;
        t_end = k * T_fundamental + 0.2;
        interval = t_start:1/fs:t_end - 1/fs;

        % Extracting noisy phase currents from Simulink
        i_d = interp1(out.noiseharmo.time, out.noiseharmo.signals.values(:,1), interval, "linear");
        i_q = interp1(out.noiseharmo.time, out.noiseharmo.signals.values(:,2), interval, "linear");

        % Identifying the fundamental frequency
        T = t_end - t_start;
        N_k = length(interval);
        Fs_k = (N_k-1)/T;
        f_k = (-N_k/2:N_k/2-1) * (Fs_k/N_k);

        % Fourier Transform
        I_d_fft = fft(i_d);  % FFT of phase A current
        I_q_fft = fft(i_q);  % FFT of phase B current

        % Identifying the fundamental frequency

        m = length(f_k);
        f = f_k(m/2+1:end);

        f = f(2:end);
        m = length(I_q_fft);
        positiveiqfft = I_q_fft(2 : m/2);

        [~, harmonique_idx] = max(abs(positiveiqfft));  % Adjust for your fundamental frequency

        % Extract the fundamental component for each phase
        I_q_fundamental = I_q_fft(1 + harmonique_idx);
        I_d_fundamental = I_d_fft(1 + harmonique_idx);

        % Calculate magnitude and phase (phasor form)
        Iq(k) = I_q_fundamental/(N_k);  % Normalize by number of samples
        Id(k) = I_d_fundamental/(N_k);
        Iq(k) = abs(Iq(k));

    end

    % Calculating statistics (rounded to 2 decimals and formatted as percentage)
    moyenne = round(mean(Iq), 4);
    ecarttype = round(std(Iq), 4);

    pd = fitdist(Iq', 'Rician');

    s = pd.s;   % Scale parameter
    sigma = pd.sigma;  % Dispersion/noise parameter

    p1 = (1 - cdf('Rician', threshold, s, sigma))*100;

    % Adding results to the table with the correct formatting
    results = [results; {Te, speed_i, Nf, sprintf('%.2f', ia), sprintf('%.2f', Rf), sprintf('%.4f', moyenne), sprintf('%.4f', ecarttype), sprintf('%.2f%%', p1)}];
end

% Creating a MATLAB table with `ia` between `Speed` and `Rf`
T = cell2table(results, 'VariableNames', {'Te[Nm]', 'Speed[rpm]', 'Nf', 'if/in', 'Rf', 'Mean[A]', 'Standard Deviation [A]', 'Detection Probability [%]'});

% Writing to an Excel file
writetable(T, 'proba.xlsx');

disp('Results have been saved in proba.xlsx');
