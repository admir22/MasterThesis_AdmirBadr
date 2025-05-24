clear; clc; close all;

% Data from the table
Te_values = [0, 0, 9.65, 9.65];
Rf = 0;
Speed_values = [500, 1000, 500, 1000];

% Initializing results
results = [];
threlistshold = zeros(1, length(Te_values));

% Loop over each case given in the table
for i = 1:length(Te_values)
    Te = Te_values(i);
    speed_i = Speed_values(i);

    % Initializing parameters
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

    unbalance_percentages = zeros(1, N_periods);

    % Calculating unbalance percentages
    for k = 1:N_periods
        t_start = (k - 1) * T_fundamental + 0.2;
        t_end = k * T_fundamental + 0.2;
        interval = t_start:1/fs:t_end - 1/fs;

        % Extracting noisy phase currents from Simulink
        i_a = interp1(out.simoutnoise.time, out.simoutnoise.signals.values(:,1), interval, "linear");
        i_b = interp1(out.simoutnoise.time, out.simoutnoise.signals.values(:,2), interval, "linear");
        i_c = interp1(out.simoutnoise.time, out.simoutnoise.signals.values(:,3), interval, "linear");

        % Fourier Transform
        I_a_fft = fft(i_a);
        I_b_fft = fft(i_b);
        I_c_fft = fft(i_c);

        % Identifying the fundamental frequency
        T = t_end - t_start;
        N_k = length(interval);
        Fs_k = (N_k-1)/T;
        f_k = (-N_k/2:N_k/2-1) * (Fs_k/N_k);
        f_k = f_k(length(f_k)/2+1:end);

        [~, fundamental_idx] = min(abs(f_k - (speed_i * 3 / 60)));

        % Extracting the fundamental components
        I_a_fundamental = I_a_fft(fundamental_idx);
        I_b_fundamental = I_b_fft(fundamental_idx);
        I_c_fundamental = I_c_fft(fundamental_idx);

        % Converting to complex numbers
        Ia = I_a_fundamental / (N_k);
        Ib = I_b_fundamental / (N_k);
        Ic = I_c_fundamental / (N_k);

        % Transforming into symmetrical components
        a = exp(1j * 2 * pi / 3);
        I0 = (1/3) * (Ia + Ib + Ic);
        I1 = (1/3) * (Ia + a*Ib + a^2*Ic);
        I2 = (1/3) * (Ia + a^2*Ib + a*Ic);

        % Calculating magnitudes and unbalance percentage
        I2_mag = abs(I2);
        I1_mag = abs(I1);
        unbalance_percentages(k) = (I2_mag / I1_mag) * 100;
    end

    % Calculating statistics (rounded to 2 decimals and formatted as percentage)
    moyenne = round(mean(unbalance_percentages), 4);
    ecarttype = round(std(unbalance_percentages), 4);

    pd = fitdist(unbalance_percentages', 'Rician');

    s = pd.s;   % Scale parameter
    sigma = pd.sigma;  % Dispersion/noise parameter

    treeshold = icdf('Rician', 0.999, s, sigma); 
    threlistshold(i) = treeshold;

    % Adding results to the table with the correct formatting
    results = [results; {Te, speed_i, sprintf('%.2f', moyenne), sprintf('%.2f', ecarttype), sprintf('%.4f', treeshold)}];
end

% Creating a MATLAB table with `ia` between `Speed` and `Rf`
T = cell2table(results, 'VariableNames', {'Te[Nm]', 'Speed[rpm]', 'Mean[%]', 'Standard Deviation [%]', 'Threshold [%]'});

% Writing to an Excel file
writetable(T, 'thresholdnegseg_closedloop.xlsx');

disp('Results have been saved in thresholdnegseg_closedloop.xlsx');
