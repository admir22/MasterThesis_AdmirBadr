% Define the values for the variables
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
If_values = [0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5, ...
             0.5, 2, 5, 0.5, 2, 5]; 

thresholdlist = [486.2434, 486.2434, 486.2434, 798.3670, 798.3670, 798.3670, ...
                 486.2434, 486.2434, 486.2434, 798.3670, 798.3670, 798.3670, ...
                 486.2434, 486.2434, 486.2434, 798.3670, 798.3670, 798.3670, ...
                 0.3657, 0.3657, 0.3657, 0.5359, 0.5359, 0.5359, ...
                 0.3657, 0.3657, 0.3657, 0.5359, 0.5359, 0.5359, ...
                 0.3657, 0.3657, 0.3657, 0.5359, 0.5359, 0.5359];

suffix_values = 1:5; % Files are numbered from 001 to 005

mu = 0;
faulty_PMSM_load;

% Initialize the results
results = [];

% Loop through all combinations of speed_i, Te, and suffix

for i = 1:length(Te_values)
    Te = Te_values(i);
    speed_i = Speed_values(i);

    t_offset = 0; % Initialize the time offset
    time = 0:1/fs:10-1/fs;
    T_fundamental = 60 / (speed_i * 3); % Fundamental period
    N_periods = floor(10 / T_fundamental);
    Iq_tot = [];
    i_q_plot = [];
    t_Iq_tot = [];
    theta_plot = [];
    ia_plot = [];
    ib_plot = [];
    ic_plot = [];

    for suffix = suffix_values
        % Dynamically construct the filename, keeping Te with a dot
        filename = sprintf('HealthyOL%d_%0.2f%03d.mat', speed_i, Te, suffix);
        data = load(filename);
        Te_str = strrep(sprintf('%.2f', Te), '.', '_'); % "9.65" → "9_65"
        filename = sprintf('HealthyOL%d_%s%03d.mat', speed_i, Te_str, suffix);

        % Remove the .mat extension to get the filename without extension
        filename_no_ext = erase(filename, '.mat');

        petit_iq = zeros(1, N_periods);
        % Calculate the imbalance percentages

        for k = 1:N_periods
            t_start = (k - 1) * T_fundamental;
            t_end = k * T_fundamental;
            interval = t_start:1/fs:t_end - 1/fs;

            % Extract noisy phase currents from Simulink
            i_q = interp1(data.(filename_no_ext).X(1).Data(:,:), data.(filename_no_ext).Y(1).Data(:,:), interval, "linear");
            thetam = thetam*2*pi/360;
            thetae = mod(thetam*p,2*pi);
            theta_plot = [theta_plot, thetae]; 
            ia_plot = [ia_plot, i_a];
            % Store all values of i_q and the corresponding time   
            i_q_plot = [i_q_plot, i_q];
            t_Iq_tot = [t_Iq_tot, interval + t_offset];

            % Identify the fundamental frequency
            T = t_end - t_start;
            N_k = length(interval);
            Fs_k = (N_k-1)/T;
            f_k = (-N_k/2:N_k/2-1) * (Fs_k/N_k);

            % Fourier Transform
            I_q_fft = fft(i_q);  % FFT of phase B current
            I_a_fft = fft(i_a);

            % Identify the fundamental frequency
            m = length(f_k);
            f = f_k(m/2+1:end);
            f = f(2:end);
            m = length(I_q_fft);
            positiveiqfft = I_q_fft(2 : m/2);

            [~, harmonique_idx] = max(abs(positiveiqfft));  % Adjust for your fundamental frequency

            % Extract the fundamental component for each phase
            I_q_fundamental = I_q_fft(1 + harmonique_idx);

            % Calculate magnitude and phase (phasor form)
            Iq = I_q_fundamental/(N_k);  % Normalize by number of samples
            petit_iq(k) = abs(Iq);
        end

        Iq_tot = [Iq_tot, petit_iq];
        Iq_tot = Iq_tot(~isnan(Iq_tot));
        t_offset = t_Iq_tot(end);
    end

    % Calculate statistics (rounded to 2 decimal places and formatted as percentages)
    moyenne = round(mean(Iq_tot), 4);
    ecarttype = round(std(Iq_tot), 4);

    pd = fitdist(Iq_tot', 'Rician');

    s = pd.s;   % Scale parameter
    sigma = pd.sigma;  % Dispersion/noise parameter

    p1 = (1 - cdf('Rician', threshold, s, sigma)) * 100; % probability

    % Add the results to the table with proper formatting
    results = [results; {Te, speed_i, Nf, sprintf('%.2f', If), sprintf('%.4f%%', moyenne), sprintf('%.4f%%', ecarttype), sprintf('%.4f%%', p1)}];

end

% Create a MATLAB table with `ia` between `Speed` and `Rf`
T = cell2table(results, 'VariableNames', {'Te[Nm]', 'Speed[rpm]', 'Nf', 'if/in', 'Mean[%]', 'Standard Deviation [%]', 'Detection Probability [%]'});

% Write to an Excel file
writetable(T, 'EXP_FAULTYclosedloopHARMO.xlsx');

disp('The results have been saved to EXP_FAULTYclosedloopHARMO.xlsx');

% Plot i_q as a function of total time
figure;
plot(t_Iq_tot, i_q_plot, 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Current i_q [A]');
title('Current i_q as a function of time');
grid on;

% Plot theta_elec as a function of total time
figure;
plot(t_Iq_tot, theta_plot, 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Theta');
title('Theta as a function of time');
grid on;

% Plot i_a as a function of total time
figure;
plot(t_Iq_tot, ia_plot, 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Current i_a [A]');
title('Current i_a as a function of time');
grid on;
