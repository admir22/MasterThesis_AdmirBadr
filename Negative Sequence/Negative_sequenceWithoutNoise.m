% Load Simulink data (assuming it's saved in a structure 'simout' with fields i_abc1, i_abc2, i_abc3)
% Replace 'simout' with the actual name of the variable if it's different.

%time = tf/2:1/fs:tf-1/fs;
time = 0:1/fs:tf-1/fs;
i_a = interp1(out.simout.time, out.simout.signals.values(:,1),time, "linear");
i_b = interp1(out.simout.time, out.simout.signals.values(:,2),time, "linear");
i_c = interp1(out.simout.time, out.simout.signals.values(:,3),time, "linear");





% Define simulation parameters
T = time(end) - time(1);  % Sampling period



% Fourier Transform
N = length(time);  % Number of samples
Fs = (N-1)/T;
f = (-N/2:N/2-1)*(Fs/N); % Frequency vector
fh = f;


I_a_fft = fft(i_a);  % FFT of phase A current
I_b_fft = fft(i_b);  % FFT of phase B current
I_c_fft = fft(i_c);  % FFT of phase C current



plot(fh,abs(fftshift(I_a_fft))/N,"LineWidth",3)
title("fft Spectrum in the Positive and Negative Frequencies")
xlabel("f (Hz)")
ylabel("|fft(X)|")
hold on


m = length(f);
f = f(m/2+1:end);



% Find fundamental frequency component
[~, fundamental_idx] = min(abs(f - (speed_i*3/60)));  % Adjust for your fundamental frequency





% Extract the fundamental component for each phase
I_a_fundamental = I_a_fft(fundamental_idx);
I_b_fundamental = I_b_fft(fundamental_idx);
I_c_fundamental = I_c_fft(fundamental_idx);

% Calculate magnitude and phase (phasor form)
Ia = I_a_fundamental/(N/2);  % Normalize by number of samples
Ib = I_b_fundamental/(N/2);
Ic = I_c_fundamental/(N/2);

% Display the results
fprintf('Phasor for I_a: Magnitude = %.3f, Angle = %.2f degrees\n', abs(Ia), angle(Ia)*180/pi);
fprintf('Phasor for I_b: Magnitude = %.3f, Angle = %.2f degrees\n', abs(Ib), angle(Ib)*180/pi);
fprintf('Phasor for I_c: Magnitude = %.3f, Angle = %.2f degrees\n', abs(Ic), angle(Ic)*180/pi);


% Transformation en composantes symétriques
a = exp(1j*2*pi/3);   % Opérateur 'a'

% Calcul des composantes de séquence positive (I1), négative (I2), et homopolaire (I0)
I0 = (1/3) * (Ia + Ib + Ic);
I1 = (1/3) * (Ia + a*Ib + a^(2)*Ic);
I2 = (1/3) * (Ia + a^2*Ib + a*Ic);

% Calcul des magnitudes
I0_mag = abs(I0);
I1_mag = abs(I1);
I2_mag = abs(I2);


% Données
variables = {'i0', 'i1', 'i2'};
valeurs = [I0_mag, I1_mag, I2_mag];

% Créer une nouvelle figure pour le graphique
figure;

% Générer le graphique en bâtonnets
bar(valeurs);

% Ajouter les étiquettes pour les axes et le titre
set(gca, 'XTickLabel', variables); % Associer les noms des variables sur l'axe x
xlabel('Variables');
ylabel('Valeurs');
title('Graphique en bâtonnets des variables ia, ib et ic');

% Ajuster l'axe y pour bien voir toutes les valeurs
ylim([0 max(valeurs)*1.2]); % Placer la limite supérieure à 1 pour bien inclure toutes les valeurs



% Vous pouvez également afficher le pourcentage de déséquilibre
unbalance_percentage = (I2_mag ./ I1_mag) * 100;
fprintf('Pourcentage de déséquilibre moyen : %f %%\n', mean(unbalance_percentage));

% Si le pourcentage de déséquilibre dépasse une certaine valeur, considérer cela comme un défaut
unbalance_threshold = 5; % Par exemple, 5%
if mean(unbalance_percentage) > unbalance_threshold
    disp('Défaut détecté : pourcentage de déséquilibre supérieur au seuil.');
else
    disp('Aucun défaut détecté selon le pourcentage de déséquilibre.');
end






% Créer des FFT remises à zéro, sauf pour la composante fondamentale
I_a_fft_fundamental = zeros(size(I_a_fft));
I_b_fft_fundamental = zeros(size(I_b_fft));
I_c_fft_fundamental = zeros(size(I_c_fft));

% Garder la composante fondamentale avec son amplitude complexe
I_a_fft_fundamental(fundamental_idx) = I_a_fft(fundamental_idx);
I_b_fft_fundamental(fundamental_idx) = I_b_fft(fundamental_idx);
I_c_fft_fundamental(fundamental_idx) = I_c_fft(fundamental_idx);

% Inverse Fourier Transform to reconstruct only the fundamental component
i_a_fundamental_reconstructed = ifft(I_a_fft_fundamental, 'symmetric');
i_b_fundamental_reconstructed = ifft(I_b_fft_fundamental, 'symmetric');
i_c_fundamental_reconstructed = ifft(I_c_fft_fundamental, 'symmetric');

% Comparer les courants originaux et reconstruits pour vérifier l'amplitude
figure;
subplot(3, 1, 1);
plot(time, i_a, 'b', time, i_a_fundamental_reconstructed, 'r--'); % Comparaison phase A
title('Phase A - Original vs Fundamental Component Only');
xlabel('Time (s)');
ylabel('Current (A)');
legend('Original', 'Reconstructed Fundamental');

subplot(3, 1, 2);
plot(time, i_b, 'b', time, i_b_fundamental_reconstructed, 'r--'); % Comparaison phase B
title('Phase B - Original vs Fundamental Component Only');
xlabel('Time (s)');
ylabel('Current (A)');
legend('Original', 'Reconstructed Fundamental');

subplot(3, 1, 3);
plot(time, i_c, 'b', time, i_c_fundamental_reconstructed, 'r--'); % Comparaison phase C
title('Phase C - Original vs Fundamental Component Only');
xlabel('Time (s)');
ylabel('Current (A)');
legend('Original', 'Reconstructed Fundamental');

% Calcul de l'erreur moyenne quadratique pour chaque phase
rmse_a = sqrt(mean((i_a - i_a_fundamental_reconstructed).^2));
rmse_b = sqrt(mean((i_b - i_b_fundamental_reconstructed).^2));
rmse_c = sqrt(mean((i_c - i_c_fundamental_reconstructed).^2));

% Affichage des erreurs
fprintf('Erreur quadratique moyenne pour la phase A : %.5f\n', rmse_a);
fprintf('Erreur quadratique moyenne pour la phase B : %.5f\n', rmse_b);
fprintf('Erreur quadratique moyenne pour la phase C : %.5f\n', rmse_c);

subplot(2,1,1);
plot(fh, abs(fftshift(I_a_fft)) / N*2, "LineWidth", 4)
title("id fft Spectrum in the Positive and Negative Frequencies", 'FontSize', 20)
xlabel("f (Hz)", 'FontSize', 20)
ylabel("|fft(X)|", 'FontSize', 20)
set(gca, 'FontSize', 20)
hold on