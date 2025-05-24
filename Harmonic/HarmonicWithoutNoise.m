


%Data from simulation
time = tf/2:1/fs:tf-1/fs;
i_d = interp1(out.simout1.time, out.simout1.signals.values(:,1),time, "linear");
i_q = interp1(out.simout1.time, out.simout1.signals.values(:,2),time, "linear");



% Define simulation parameters
T = time(end) - time(1);  % Sampling period



% Fourier Transform
N = length(time);  % Number of samples
Fs = (N-1)/T;
fh = (-N/2:N/2-1)*(Fs/N); % Frequency vector



I_d_fft = fft(i_d);  % FFT of d-component of the current
I_q_fft = fft(i_q);  % FFT of q-component of the current

I_d_fft(1) = 0;
I_q_fft(1) = 0;




m = length(fh);
f = fh(m/2+1:end);

f = f(2:end); %Remove DC Component


m = length(I_q_fft);
positiveiqfft = I_q_fft(2 : m/2);


[~, harmonique_idx] = max(abs(positiveiqfft)); %Take biggest harmonic index



% Extract the harmonic component for d,q currents
I_q_fundamental = I_q_fft(1 + harmonique_idx);
I_d_fundamental = I_d_fft(1 + harmonique_idx);

% Calculate magnitude and phase (phasor form)
Iq = I_q_fundamental/(N);  % Normalize by number of samples/2
Id = I_d_fundamental/(N);

if abs(Iq) > 0.5
    disp('Fault detected: Harmonic above threshold.');
else
    disp('No fault detected.');
end

fprintf('Amplitude of the biggest iq harmonics : %f A\n', abs(Iq));

fprintf('Amplitude of the 2biggest iq harmonics : %f A\n', abs(Iq));


% ----- Plot FFT for i_d and i_q -----
figure;
set(gcf, 'Color', 'w'); % Fond blanc

% Options graphiques
line_width = 5;
font_size = 35;  % Taille de la police
legend_size = 30;
axis_thickness = 5;

% Plot d-axis current (id) FFT
subplot(2,1,1);
plot(fh, abs(fftshift(I_d_fft)) / N, 'LineWidth', line_width);
xlim([50, 250]); % Plage de fréquence de [-200, 200] Hz
title('id FFT Spectrum', 'FontSize', font_size);
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('|FFT(i_d)|', 'FontSize', font_size);
set(gca, 'FontSize', font_size, 'LineWidth', axis_thickness); % Améliorer la visibilité des axes
grid on;

% Plot q-axis current (iq) FFT
subplot(2,1,2);
plot(fh, abs(fftshift(I_q_fft)) / N, 'LineWidth', line_width);
xlim([50, 250]); % Plage de fréquence de [-200, 200] Hz
title('iq FFT Spectrum', 'FontSize', font_size);
xlabel('Frequency (Hz)', 'FontSize', font_size);
ylabel('|FFT(i_q)|', 'FontSize', font_size);
set(gca, 'FontSize', font_size, 'LineWidth', axis_thickness); % Améliorer la visibilité des axes
grid on;

