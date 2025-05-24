% Define the values for the variables
Te_values = [0, 0, 9.65, 9.65]; % List of Te values
Speed_values = [500, 1000, 500, 1000]; % List of speed_i values

mu = 0;
faulty_PMSM_load;

% Initialize the results
results = [];
threlistshold = zeros(1, length(Te_values));

% Loop through all combinations of speed_i, Te, and suffix

for i = 1:length(Te_values)
    Te = Te_values(i);
    speed_i = Speed_values(i);

    if speed_i == 500
        suffix_values = 1:10; % Files are numbered from 001 to 005
    else
        suffix_values = 1:5; % Files are numbered from 001 to 005
    end

    time = 0:1/fs:10-1/fs;
    T_fundamental = 60 / (speed_i * 3); % Fundamental period
    N_periods = floor(10 / T_fundamental);
    unbalance_percentages = [];
    ia_bruit = [];

    for suffix = suffix_values
        % Dynamically construct the filename, keeping Te with a dot
        filename = sprintf('HealthyOL%d_%0.2f%03d.mat', speed_i, Te, suffix);
        data = load(filename);
        Te_str = strrep(sprintf('%.2f', Te), '.', '_'); % "9.65" → "9_65"
        filename = sprintf('HealthyOL%d_%s%03d.mat', speed_i, Te_str, suffix);

        % Remove the .mat extension to get the filename without extension
        filename_no_ext = erase(filename, '.mat');

        petit_unb_per = zeros(1, N_periods);
        % Calculate the imbalance percentages

        for k = 1:N_periods
            t_start = (k - 1) * T_fundamental;
            t_end = k * T_fundamental;
            interval = t_start:1/fs:t_end - 1/fs;

            % Extract noisy phase currents from Simulink
            i_a = interp1(data.(filename_no_ext).X(1).Data(:,:), data.(filename_no_ext).Y(1).Data(:,:), interval, "linear");
            i_b = interp1(data.(filename_no_ext).X(1).Data(:,:), data.(filename_no_ext).Y(2).Data(:,:), interval, "linear");
            i_c = interp1(data.(filename_no_ext).X(1).Data(:,:), data.(filename_no_ext).Y(3).Data(:,:), interval, "linear");

            % Fourier Transform
            I_a_fft = fft(i_a);
            I_b_fft = fft(i_b);
            I_c_fft = fft(i_c);
            ia_bruit = [ia_bruit, i_a];
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
            petit_unb_per(k) = I2_mag;
        end

        unbalance_percentages = [unbalance_percentages, petit_unb_per];

        unbalance_percentages = unbalance_percentages(~isnan(unbalance_percentages));

    end

    % Calculate statistics (rounded to 2 decimal places and formatted as percentages)
    moyenne = round(mean(unbalance_percentages), 2);
    ecarttype = round(std(unbalance_percentages), 2);

    figure

    pd = fitdist(unbalance_percentages', 'Rician');

    histogram(unbalance_percentages', 'Normalization', 'pdf', 'FaceAlpha', 0.3, 'FaceColor', 'r');
    hold on

    s = pd.s;   % Scale parameter
    sigma = pd.sigma;  % Dispersion/noise parameter

    if Te < 1
        % Define x values for plotting
        x = 0:0.00001:0.02; % Range of values
    else
        x = 0:0.0001:0.2;
    end

    % Calculate the probability density function of the Rician distribution
    y = pdf('Rician', x, s, sigma);

    plot(x, y, 'r', 'LineWidth', 3);

    hold on

    % Calculate the threshold where P(X >= threshold) = 5% (i.e., P(X <= threshold) = 95%)
    treeshold = icdf('Rician', 0.999, s, sigma);

    plot([treeshold treeshold], ylim, 'm--', 'LineWidth', 2, 'DisplayName', sprintf('threshold: %.2f', treeshold));

    threlistshold(i) = treeshold;

    % Add the results to the table with proper formatting
    results = [results; {Te, speed_i, sprintf('%.2f%%', moyenne), sprintf('%.2f%%', ecarttype), sprintf('%.2f%%', treeshold)}];

end

% Create a MATLAB table with `ia` between `Speed` and `Rf`
T = cell2table(results, 'VariableNames', {'Te[Nm]', 'Speed[rpm]', 'Mean[%]', 'Standard Deviation [%]', 'Threshold [%]'});

% Write to an Excel file
writetable(T, 'EXP_thresholdOPENloopnegseq.xlsx');

disp('The results have been saved to EXP_thresholdOPENloopnegseq.xlsx');
