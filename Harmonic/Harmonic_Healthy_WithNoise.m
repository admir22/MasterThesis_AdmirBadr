clear; clc; close all;

% Data from the table
Te_values = [0, 0, 9.65, 9.65];
Rf = 0;
Speed_values = [500, 1000, 500, 1000];

% Initializing results
results = [];
threlistshold = zeros(1, length(Te_values));

% Loop over each case provided in the table
for i = 1:length(Te_values)
    Te = Te_values(i);
    speed_i = Speed_values(i);

    % Initialize parameters
    mu = (0/378); % Fault severity
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

    out = sim('faulty_motor_model_closedloop');

    % Simulation parameters
    time = 0.2:1/fs:tf-1/fs; 
    T_fundamental = 60 / (speed_i * 3); % Fundamental period
    N_periods = floor((time(end)-time(1)) / T_fundamental);

    Iq = zeros(1, N_periods);
    Id = zeros(1, N_periods);

    % Calculation of unbalance percentages
    for k = 1:N_periods
        t_start = (k - 1) * T_fundamental + 0.2;
        t_end = k * T_fundamental + 0.2;
        interval = t_start:1/fs:t_end - 1/fs;

        % Extract noisy phase currents from Simulink
        i_d = interp1(out.noiseharmo.time, out.noiseharmo.signals.values(:,1), interval, "linear");
        i_q = interp1(out.noiseharmo.time, out.noiseharmo.signals.values(:,2), interval, "linear");

        % Identification of fundamental frequency
        T = t_end - t_start;
        N_k = length(interval);
        Fs_k = (N_k-1)/T;
        f_k = (-N_k/2:N_k/2-1) * (Fs_k/N_k);

        % Fourier Transform
        I_d_fft = fft(i_d);  % FFT of phase A current
        I_q_fft = fft(i_q);  % FFT of phase B current

        % Identification of the fundamental frequency

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

    % Calculating statistics (rounded to 4 decimal places and formatted as percentages)
    moyenne = round(mean(Iq), 4);
    ecarttype = round(std(Iq), 4);

    pd = fitdist(Iq', 'Rician');

    s = pd.s;   % Scale parameter
    sigma = pd.sigma;  % Dispersion/noise parameter

    treeshold = icdf('Rician', 0.999, s, sigma); 
    threlistshold(i) = treeshold;

    % Adding results to the table with the correct formatting
    results = [results; {Te, speed_i, sprintf('%.4f', moyenne), sprintf('%.4f', ecarttype), sprintf('%.4f', treeshold)}];
end

% Creating a MATLAB table with `ia` between `Speed` and `Rf`
T = cell2table(results, 'VariableNames', {'Te[Nm]', 'Speed[rpm]', 'Mean[A]', 'Standard Deviation [A]', 'Threshold [A]'});

% Writing to an Excel file
writetable(T, 'thresholdclosedloopharmonics.xlsx');

disp('Results have been saved in thresholdclosedloopharmonics.xlsx');
