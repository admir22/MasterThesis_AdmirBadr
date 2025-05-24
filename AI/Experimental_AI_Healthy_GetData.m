% Define the values for the variables
Te_values = [0, 0, 0, 0, 0, 0, ...
             0, 0, 0, 0, 0, 0, ...
             0, 0, 0, 0, 0, 0, ...
             9.65, 9.65, 9.65, 9.65, 9.65, 9.65, ...
             9.65, 9.65, 9.65, 9.65, 9.65, 9.65, ...
             9.65, 9.65, 9.65, 9.65, 9.65, 9.65];

Nf_values = [5, 10, 10, 1, 4, 10, ...
             5, 20, 25, 5, 5, 25, ...
             10, 20, 25, 10, 10, 25, ...
             2, 10, 10, 1, 4, 10, ...
             5, 20, 25, 5, 5, 25, ...
             10, 20, 25, 10, 10, 25];

If_values = [0.1, 0.2, 0.5, 0.5, 0.5, 0.5, ...
             0.5, 0.5, 1, 0.5, 2, 1, ...
             0.5, 1, 2, 0.5, 2, 2, ...
             0.25, 0.2, 0.5, 0.5, 0.5, 0.5, ...
             0.5, 0.5, 1, 0.5, 2, 1, ...
             0.5, 1, 2, 0.5, 2, 2]; 

Speed_values = [500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000, ...
                500, 500, 500, 1000, 1000, 1000];

mu = 0;
faulty_PMSM_load;

% Initialize the results
results = [];

% Simulation parameters
Te_total = []; % To store all Te values
Nf_total = []; % To store all Nf values
Speed_total = []; % To store all speeds
Ia_total = []; % To store all Ia values
Iq_2nd_inst_total = []; % To store all 2nd harmonic values
Iq_4th_inst_total = []; % To store all 4th harmonic values

% Loop through all combinations of Te, Speed_i, and suffix
for i = 1:length(Te_values)
    Te = Te_values(i);
    speed_i = Speed_values(i);
    If = If_values(i);
    Nf = Nf_values(i);

    if speed_i == 1000
        suffix_values = 1:5; % Files are numbered from 001 to 005
    else
        suffix_values = 1:10; % Files are numbered from 001 to 005
    end

    % Construct the overall time series for 50 seconds
    t_offset = 0; % Initialize the time offset
    T_fundamental = 60 / (speed_i * 3); % Fundamental period
    N_periods = floor((10) / T_fundamental); % Period over 10 seconds

    % Loop over each file (suffix from 001 to 005)
    for suffix = suffix_values
        % Dynamically construct the filename, keeping Te with a dot
        filename = sprintf('OLFaulty%d_%.2f_%d_%.2f%03d.mat', speed_i, Te, Nf, If, suffix);
        data = load(filename);
        Te_str = strrep(sprintf('%.2f', Te), '.', '_'); % "9.65" → "9_65"
        If_str = strrep(sprintf('%.2f', If), '.', '_'); % "9.65" → "9_65"
        filename = sprintf('OLFaulty%d_%s_%d_%s%03d.mat', speed_i, Te_str, Nf, If_str, suffix);

        % Remove the .mat extension to get the filename without extension
        filename_no_ext = erase(filename, '.mat');

        % Extract data from the file
        time = 0.2:1/fs:tf-1/fs;
        T_fundamental = 60 / (speed_i * 3);
        N_periods = floor((time(end)-time(1)) / T_fundamental);

        % Loop over each period to calculate the harmonics
        for k = 1:N_periods
            t_start = (k - 1) * T_fundamental;
            t_end = k * T_fundamental;
            interval = t_start:1/fs:t_end - 1/fs;

            % Extract the noisy phase currents
            i_q = interp1(data.(filename_no_ext).X(1).Data(:,:), data.(filename_no_ext).Y(6).Data(:,:), interval, "linear");

            N = length(interval);

            fh = (-N/2:N/2-1)*(Fs/N); % Frequency vector

            % Calculate harmonics
            I_q_fft = fft(i_q);
            positiveiqfft = I_q_fft(2:end/2);
            [~, harmonique_idx] = max(abs(positiveiqfft));
            harmonique_idx4 = 2 * harmonique_idx;

            Iq_2nd_inst = 20 * log10(abs(I_q_fft(1 + harmonique_idx)) / length(interval));
            Iq_4th_inst = 20 * log10(abs(I_q_fft(1 + harmonique_idx4)) / length(interval));

            % Concatenate results over total time
            Te_total = [Te_total, repmat(Te, 1, length(interval))];
            Nf_total = [Nf_total, repmat(Nf, 1, length(interval))];
            Speed_total = [Speed_total, repmat(speed_i, 1, length(interval))];
            Ia_total = [Ia_total, repmat(If, 1, length(interval))]; % Ensure Ia is defined
            Iq_2nd_inst_total = [Iq_2nd_inst_total, repmat(Iq_2nd_inst, 1, length(interval))];
            Iq_4th_inst_total = [Iq_4th_inst_total, repmat(Iq_4th_inst, 1, length(interval))];
        end
    end
end

% Group all data into a final table
results = [Te_total', Nf_total', Speed_total', Ia_total', Iq_2nd_inst_total', Iq_4th_inst_total'];

% Create a MATLAB table with the results
T = table(Te_total', Nf_total', Speed_total', Ia_total', Iq_2nd_inst_total', Iq_4th_inst_total', ...
    'VariableNames', {'Te', 'Nf', 'Speed', 'Ia', 'Iq_2nd_inst', 'Iq_4th_inst'});

%AD = table2array(T);

% Save the results
%save('harmonics_data_instantaneousEXPAD.mat', 'AD');

%disp('The results have been saved in merged_results.xlsx');
