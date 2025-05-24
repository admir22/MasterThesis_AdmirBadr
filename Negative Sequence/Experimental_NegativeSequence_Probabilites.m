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

suffix_values = 1:5; % The files are numbered from 001 to 005

mu=0;
faulty_PMSM_load;

% Initialize the results
results = [];

% Loop through all combinations of speed_i, Te, and suffix

for i = 1:length(Te_values)
    Te = Te_values(i);
    speed_i = Speed_values(i);
    Nf = Nf_values(i);
    If = If_values(i);
    threshold = thresholdlist(i);

    time = 0:1/fs:10-1/fs;
    T_fundamental = 60 / (speed_i * 3); % Fundamental period
    N_periods = floor(10 / T_fundamental);
    unbalance_percentages = [];

    for suffix = suffix_values
            % Construire dynamiquement le nom du fichier en gardant Te avec un point
        filename_faulty = sprintf('OLFaulty%d_%.2f_%d_%.2f%03d.mat', speed_i, Te, Nf, If, suffix);
        data_f = load(filename_faulty);
        Te_str = strrep(sprintf('%.2f', Te), '.', '_'); % "9.65" → "9_65"
        If_str = strrep(sprintf('%.2f', If), '.', '_'); % "9.65" → "9_65"
        filename_faulty = sprintf('OLFaulty%d_%s_%d_%s%03d.mat', speed_i, Te_str, Nf, If_str, suffix);

        name_f = erase(filename_faulty, '.mat');
        if length(data_f.(name_f).Y) < 15, continue; end

        % Remove the .mat extension to get the filename without extension
        filename_no_ext = erase(filename, '.mat');

        petit_unb_per = zeros(1, N_periods);
        
        % Calculate the imbalance percentages
        for k = 1:N_periods
            t_start = (k - 1) * T_fundamental;
            t_end = k * T_fundamental;
            interval = t_start:1/fs:t_end - 1/fs;

            % Extract the noisy phase currents from Simulink
            i_a = interp1(data.(filename_no_ext).X(1).Data(:,:), data.(filename_no_ext).Y(1).Data(:,:), interval, "linear");
            i_b = interp1(data.(filename_no_ext).X(1).Data(:,:), data.(filename_no_ext).Y(2).Data(:,:), interval, "linear");
            i_c = interp1(data.(filename_no_ext).X(1).Data(:,:), data.(filename_no_ext).Y(3).Data(:,:), interval, "linear");

            % Fourier Transform
            I_a_fft = fft(i_a);
            I_b_fft = fft(i_b);
            I_c_fft = fft(i_c);

            % Identify the fundamental frequency
            T = t_end - t_start;
            N_k = length(interval);
            Fs_k = (N_k-1)/T;
            f_k = (-N_k/2:N_k/2-1) * (Fs_k/N_k);
            f_k = f_k(length(f_k)/2+1:end);

            [~, fundamental_idx] = min(abs(f_k - (speed_i * 3 / 60)));

            % Extract the fundamental components
            I_a_fundamental = I_a_fft(fundamental_idx);
            I_b_fundamental = I_b_fft(fundamental_idx);
            I_c_fundamental = I_c_fft(fundamental_idx);

            % Convert to complex numbers
            Ia = I_a_fundamental / (N_k);
            Ib = I_b_fundamental / (N_k);
            Ic = I_c_fundamental / (N_k);

            % Transform to symmetrical components
            a = exp(1j * 2 * pi / 3);
            I0 = (1/3) * (Ia + Ib + Ic);
            I1 = (1/3) * (Ia + a*Ib + a^2*Ic);
            I2 = (1/3) * (Ia + a^2*Ib + a*Ic);

            % Calculate magnitudes and imbalance percentage
            I2_mag = abs(I2);
            I1_mag = abs(I1);
            petit_unb_per(k) = (I2_mag / I1_mag) * 100;
        end

        unbalance_percentages = [unbalance_percentages, petit_unb_per];
        unbalance_percentages = unbalance_percentages(~isnan(unbalance_percentages));
    end

    % Calculate statistics (rounded to 2 decimal places and formatted as percentages)
    moyenne = round(mean(unbalance_percentages), 2);
    ecarttype = round(std(unbalance_percentages), 2);

    pd = fitdist(unbalance_percentages', 'Rician');
    s = pd.s;   % Scale parameter
    sigma = pd.sigma;  % Dispersion/noise parameter

    p1 = (1 - cdf('Rician', threshold, s, sigma))*100; % probability

    % Add the results to the table with proper formatting
    results = [results; {Te, speed_i, Nf, sprintf('%.2f', If), sprintf('%.2f%%', moyenne), sprintf('%.2f%%', ecarttype), sprintf('%.2f%%', p1)}];

end

% Create a MATLAB table with `ia` between `Speed` and `Rf`
T = cell2table(results, 'VariableNames', {'Te[Nm]', 'Speed[rpm]', 'Nf', 'if/in', 'Mean[%]', 'Standard Deviation [%]', 'Detection Probability [%]'});

% Write to an Excel file
writetable(T, 'EXP_FAULTYclosedloopnegseq.xlsx');

disp('The results have been saved to EXP_FAULTYclosedloopnegseq.xlsx');
