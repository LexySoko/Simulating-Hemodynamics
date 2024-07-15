% Fetch rawF data for violet light
session_key.subject_fullname = '...';
session_key.session_date = '...';
session_key.session_number = ...;
param_set_ids = getPrimaryWFPreprocParams;
manual_mask = fetch1(widefield.ManualMask & session_key & param_set_ids, 'mask_binned');
manual_mask = reshape(manual_mask, [], 1);
num_pixel = 10;   %adjust if needed
pixel_selected = randsample(find(manual_mask), num_pixel);
pixel_filter = gen_filter('Pixel_id', pixel_selected);
rawf = fetchn(widefield.RawF & session_key & param_set_ids & pixel_filter, 'rawf');
rawf = cell2mat(rawf)';
strobe_sequence = fetch1(widefield.StrobeSeq & session_key & param_set_ids, 'strobe_sequence');
rawf_violet = rawf(strobe_sequence == 2, :);

% Compute mode-based baseline for violet
f0_violet_mode = computeF0mode(rawf_violet, 30, 10, false);

% Compute fractional fluorescence for violet
fractional_f_violet = rawf_violet ./ f0_violet_mode;

% Pathlength for violet light
pathlength_violet = 0.022148; %mm  
molar_extinction_violet = 244764; % cm^-1 M^-1

%Above values from tabulated supplementary materials in “Wide-field optical mapping of neural activity and brain haemodynamics: considerations and novel approaches”

.......................................................................................................................................................................................

% Initialize absorption change matrix
absorption_change = zeros(size(fractional_f_violet));
change_hbt = zeros(size(fractional_f_violet));

% Iterate through selected pixels and time to compute absorption change using equation 2.6 from "Wide-field optical mapping of neural activity and brain haemodynamics:considerations and novel approaches
for pixel = 1:num_pixel
    for t = 1:size(fractional_f_violet, 1)
        f_f0 = fractional_f_violet(t, pixel);  % Fetch fractional fluorescence for pixel at time t
        log_fluorescence = log(f_f0);  % Take natural log of (f/f0)
        x = -1 / pathlength_violet;  % Calculate x
        absorption_change(t, pixel) = x * log_fluorescence;  % Compute absorption change
        change_hbt(t, pixel) = absorption_change(t, pixel) / molar_extinction_violet;  % Compute change in total hemoglobin
    end
end

% Plot change in total hemoglobin and raw fluorescence for one pixel to verify change in total hemoglobin was calculated correctly. Plots should be inverse of each other. 
figure;
pixel = 1;  % Select the first pixel for plotting

% Shorten the time window
time_window = 1:1000;  % Adjust this range as needed for a clear view of spikes

% Plot raw fluorescence
subplot(2, 1, 1);
plot(time_window, rawf_violet(time_window, pixel));
title(['Raw Fluorescence for Pixel ', num2str(pixel)]);
xlabel('Time');
ylabel('Raw Fluorescence');
grid on;

% Plot change in total hemoglobin
subplot(2, 1, 2);
plot(time_window, change_hbt(time_window, pixel), 'k');
title(['Change in Total Hemoglobin (Δ[HbT]) (Violet Light) for Pixel ', num2str(pixel)]);
xlabel('Time');
ylabel('Change in HbT (M)');
grid on;

.......................................................................................................................................................................................

%Least Squares Method (Done per pixel) for finding a hemodynamic response function. Values and equations from "Resting-state hemodynamics are spatiotemporally coupled to synchronized and symmetric neural activity in excitatory neurons”
% Step 1: Calculate Delta F / F 
delta_f_over_f_violet = (rawf_violet - f0_violet_mode) ./ f0_violet_mode;

% Step 2: Form the System Input Matrix X
n = 5000; % Length of the HRF. Adjust if needed. Larger n values are more computationally burdensome.
X = zeros(size(delta_f_over_f_violet, 1), n);

for i = 1:size(delta_f_over_f_violet, 1)
    for j = 1:min(i, n)
        X(i, j) = delta_f_over_f_violet(i - j + 1, 1);  % Using first pixel as example
    end
end

% Step 3: Form the System Output Vector y
change_hbt_pixel = change_hbt(:, 1);  % Using first pixel as example

% Step 4: Define the lambda regularization parameter
lambda = 0.01;  

% Step 5: Perform Least-Square Deconvolution to find the HRF
HRF = (X' * X + lambda * eye(n)) \ (X' * change_hbt_pixel); 

% Step 6: Plot the Estimated HRF
figure;
plot(HRF);
title('Estimated Hemodynamic Response Function (HRF)');
xlabel('Time');
ylabel('HRF Amplitude');
grid on;

% Step 7: Calculate Pearson's Correlation Coefficient
y_predicted = X * HRF;
correlation_coefficient = corr(change_hbt_pixel, y_predicted);
disp(['Pearson Correlation Coefficient: ', num2str(correlation_coefficient)]);


