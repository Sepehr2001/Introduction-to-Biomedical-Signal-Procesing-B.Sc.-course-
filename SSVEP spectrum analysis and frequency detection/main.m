clc
clear
close all
load('7Data.mat')

%% 1) plot raw data.
% Plot side by side with tiledlayout
data = X{1, 1};
Fs = 500;
stim_freq = 9.25;
% Time vector
N = length(data);
t = (0:N-1)/Fs;

% FFT
Y = fft(data);
P2 = abs(Y/N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:floor(N/2)) * (Fs/N);

figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Time domain
nexttile;
plot(t, data, 'k', 'LineWidth', 1.1);
grid on;
xlim([t(1) t(end)]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Raw Signal');

% FFT
nexttile;
plot(f, P1, 'k', 'LineWidth', 1.1);
grid on;
xlim([0 40]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title(['FFT - ' num2str(stim_freq) ' Hz']);
xline(stim_freq, 'r--', 'LineWidth', 1);
xline(2*stim_freq, 'b--', 'LineWidth', 1);

sgtitle('Raw Signal and FFT');

%% 2) find Fs

idx = 2;
data = X{idx, 1};
% Parameters
% Fs = 500;
Fs = 256;
stim_freq = 9.25;
%  bandpass filtering first
bp = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 6, ...
    'HalfPowerFrequency2', 30, ...
    'SampleRate', Fs);
data_filt = filtfilt(bp, data);
% FFT of filtered signal
N = length(data_filt);
Y = fft(data_filt);
P2 = abs(Y/N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:floor(N/2)) * (Fs/N);
% Plot
figure;
plot(f, P1, 'b-', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title(['FFT of Bandpass-Filtered Signal, Fs = ', num2str(Fs)]);
xlim([0 50]);
grid on;
xline(stim_freq, 'r--', [num2str(stim_freq), ' Hz'], 'LineWidth', 1.2);
xline(2*stim_freq, 'g--', [num2str(2*stim_freq), ' Hz'], 'LineWidth', 1.2);

%% 3) sort the 12 SSVEP signals based on the stimulus frequencies

% Manually correct the order of input SSVEP data
data_sorted = zeros(length(X{1,1}), 12);
data_sorted(:, 1) = X{1, 1};
data_sorted(:, 2) = X{4, 1};
data_sorted(:, 3) = X{7, 1};
data_sorted(:, 4) = X{10, 1};
data_sorted(:, 5) = X{2, 1};
data_sorted(:, 6) = X{5, 1};
data_sorted(:, 7) = X{8, 1};
data_sorted(:, 8) = X{11, 1};
data_sorted(:, 9) = X{3, 1};
data_sorted(:, 10) = X{6, 1};
data_sorted(:, 11) = X{9, 1};
data_sorted(:, 12) = X{12, 1};

%% 4) filter and visualize
% Parameters
Fs = 256;
stim_freq_range = 9.25:0.5:14.75;

% Bandpass filter
bp = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 6, ...
    'HalfPowerFrequency2', 30, ...
    'SampleRate', Fs);

figure;

t = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

for idx = 1:12
    data = data_sorted(:, idx);
    stim_freq = stim_freq_range(idx);
    
    % Filter signal
    data_filt = filtfilt(bp, data);
    
    % FFT
    N = length(data_filt);
    Y = fft(data_filt);
    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (0:floor(N/2)) * (Fs/N);
    
    % Tile plot
    ax = nexttile;
    h1 = plot(ax, f, P1, 'k-', 'LineWidth', 1.2);
    hold(ax, 'on');
    h2 = xline(ax, stim_freq, 'r--', 'LineWidth', 1);
    h3 = xline(ax, 2*stim_freq, 'b--', 'LineWidth', 1);
    hold(ax, 'off');
    
    xlim(ax, [0 40]);
    grid(ax, 'on');
    title(ax, ['SSVEP',num2str(idx),'~~', num2str(stim_freq) ' Hz']);
    xlabel(ax, 'f (Hz)');
    ylabel(ax, 'Magnitude');
end

title(t, 'Filterd SSVEP FFT Spectra for All 12 Stimulus Frequencies');
legend([h1 h2 h3], {'FFT','1^{st} harmonic','2^{nd} harmonic'}, ...
    'Location','northeast');

%% 5.1) AR model order choose

% Maximum AR order to test
max_order = 300;

% Create 4x3 subplot for AIC curves
figure;
t = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');

orders = 1:max_order;

for idx = 1:12
    data = data_sorted(:, idx);
    stim_freq = stim_freq_range(idx);
    
    % Filter signal
    data_filt = filtfilt(bp, data);
    
    % Compute AIC values
    aic_values = compute_aic_values_AR(data_filt, max_order);
    
    % Find optimal order
    [min_aic, opt_order] = min(aic_values);
    
    % Plot
    ax = nexttile;
    plot(ax, orders, aic_values, 'b-', 'LineWidth', 1.2);
    hold(ax, 'on');
    plot(ax, opt_order, min_aic, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    grid(ax, 'on');
    xlabel(ax, 'AR Order');
    ylabel(ax, 'AIC');
    title(ax, ['Freq: ' num2str(stim_freq) ' Hz, Optimal p=' num2str(opt_order)]);
    xlim(ax, [1 max_order]);
end

title(t, 'AIC vs AR Order for All 12 SSVEP Signals');

%% 5.2) plot the AR psd estimation using the optimal orders
% Bandpass filter
bp = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 6, ...
    'HalfPowerFrequency2', 30, ...
    'SampleRate', Fs);

% Create 3x4 subplot
figure;
t = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

for idx = 1:12
    data = data_sorted(:, idx);
    stim_freq = stim_freq_range(idx);
    
    % Filter signal
    data_filt = filtfilt(bp, data);
    
    % Find optimal AR order (testing 1 to 300)
    aic_values = compute_aic_values_AR(data_filt, 300);
    [~, opt_order] = min(aic_values);       % very high order
    
    % AR model estimation (Burg method)
    % opt_order = 10;      % very low order
    opt_order = 40;         % optimal order based on AIC curve
    [a, e] = arburg(data_filt, opt_order);
    
    % Compute AR PSD manually
    nfft = 512;
    H = freqz(1, a, nfft, Fs);
    Sxx_ar = e * abs(H).^2 / nfft;
    f_ar = (0:nfft/2)' * (Fs/nfft);
    Sxx_ar = Sxx_ar(1:nfft/2+1);
    
    % Plot
    ax = nexttile;
    h1 = plot(ax, f_ar, 10*log10(Sxx_ar), 'k-', 'LineWidth', 1.2);
    grid(ax, 'on');
    xlim(ax, [0 40]);
    xlabel(ax, 'Freq (Hz)');
    ylabel(ax, 'PSD (dB/Hz)');
    title(ax, ['AR PSD (p=' num2str(opt_order) '): ' num2str(stim_freq) ' Hz']);
    h2 = xline(ax, stim_freq, 'r--', 'LineWidth', 1);
    h3 = xline(ax, 2*stim_freq, 'b--', 'LineWidth', 1);
end

title(t, 'AR Model PSD with Optimal Order (AIC) for All 12 Frequencies');
legend([h1 h2 h3], {'FFT','1^{st} harmonic','2^{nd} harmonic'}, ...
    'Location','northeast');


%% 6.1) MA model order choose
max_order = 50;

% Filter all signals first
data_filt_all = zeros(size(data_sorted));
for idx = 1:12
    data_filt_all(:,idx) = filtfilt(bp, data_sorted(:,idx));
end

% Start process-based pool
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');
end

% Preallocate BIC matrix
bic_matrix = inf(12, max_order);

% Parallel computation
parfor idx = 1:12
    x = data_filt_all(:, idx);
    N = length(x);
    local_bic = inf(1, max_order);

    for q = 1:max_order
        try
            Mdl = arima(0,0,q);
            EstMdl = estimate(Mdl, x, 'Display', 'off');
            [~,~,logL] = infer(EstMdl, x);
            k = q + 1;   % q MA coeffs + variance
            local_bic(q) = -2*logL + k*log(N);
        catch
            local_bic(q) = inf;
        end
    end

    bic_matrix(idx,:) = local_bic;
end

figure;
t = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');

for idx = 1:12
    stim_freq = stim_freq_range(idx);
    bic_values = bic_matrix(idx,:);
    bic_values = bic_values(:).';
    orders = 1:length(bic_values);

    valid_idx = isfinite(bic_values);

    ax = nexttile;

    if any(valid_idx)
        valid_orders = orders(valid_idx);
        valid_bic = bic_values(valid_idx);

        [min_bic, local_idx] = min(valid_bic);
        opt_order = valid_orders(local_idx);

        plot(ax, valid_orders, valid_bic, 'b-', 'LineWidth', 1.2);
        hold(ax, 'on');
        plot(ax, opt_order, min_bic, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
        title(ax, ['Freq: ' num2str(stim_freq) ' Hz, q=' num2str(opt_order)]);
    else
        title(ax, ['Freq: ' num2str(stim_freq) ' Hz - no valid fit']);
    end

    grid(ax, 'on');
    xlabel(ax, 'MA Order');
    ylabel(ax, 'BIC');
    xlim(ax, [1 max_order]);
end

title(t, 'BIC vs MA Order for All 12 SSVEP Signals');

%% 6.2) Plot MA PSD estimation using the optimal orders
figure;
t = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

nfft = 512;
opt_orders = zeros(12,1);

for idx = 1:12
    data_filt = data_filt_all(:, idx);
    stim_freq = stim_freq_range(idx);

    bic_values = bic_matrix(idx,:);
    valid_idx = isfinite(bic_values);

    ax = nexttile;

    if ~any(valid_idx)
        title(ax, [num2str(stim_freq) ' Hz - no valid fit']);
        axis(ax, 'off');
        continue
    end

    valid_orders = find(valid_idx);
    valid_bic = bic_values(valid_idx);
    [~, local_idx] = min(valid_bic);
    opt_order = valid_orders(local_idx);
    opt_orders(idx) = opt_order;
    % opt_order =12;  % optimal order based on BIC curve
    try
        Mdl = arima(0,0,opt_order);
        EstMdl = estimate(Mdl, data_filt, 'Display', 'off');

        ma_coeffs = cell2mat(EstMdl.MA);
        b = [1 ma_coeffs];
        sigma2 = EstMdl.Variance;

        H = freqz(b, 1, nfft, Fs);
        Sxx_ma = sigma2 * abs(H).^2;

        f_ma = linspace(0, Fs/2, nfft/2 + 1);
        Sxx_ma = Sxx_ma(1:nfft/2+1);

        h1 = plot(ax, f_ma, 10*log10(Sxx_ma), 'k-', 'LineWidth', 1.2);
        hold(ax, 'on');
        h2 = xline(ax, stim_freq, 'r--', 'LineWidth', 1);
        h3 = xline(ax, 2*stim_freq, 'b--', 'LineWidth', 1);

        grid(ax, 'on');
        xlim(ax, [0 40]);
        xlabel(ax, 'Freq (Hz)');
        ylabel(ax, 'PSD (dB)');
        title(ax, ['MA PSD (q=' num2str(opt_order) '): ' num2str(stim_freq) ' Hz']);

    catch
        title(ax, [num2str(stim_freq) ' Hz - estimation failed']);
        axis(ax, 'off');
    end
end

title(t, 'MA Model PSD with Optimal Order (BIC) for All 12 Frequencies');
legend([h1 h2 h3], {'MA PSD','1^{st} harmonic','2^{nd} harmonic'}, ...
    'Location','northeast');

%% 7.1) ARMA model order choose

% ARMA candidate orders: [p q]
arma_orders = [
    1 1
    2 3
    4 5
    6 8
    8 9
    12 3
    16 2
    25 6
    40  10
    50  12
];

num_models = size(arma_orders,1);
aic_matrix = inf(12, num_models);
best_pq = nan(12,2);

% Optional parallel computation
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');
end

parfor idx = 1:12
    x = data_filt_all(:,idx);
    N = length(x);
    local_aic = inf(1, num_models);

    for k = 1:num_models
        p = arma_orders(k,1);
        q = arma_orders(k,2);

        try
            Mdl = arima(p,0,q);
            EstMdl = estimate(Mdl, x, 'Display', 'off');

            % Use innovation variance as E_pq
            Epq = EstMdl.Variance;

            % AIC(p,q) = N*ln(Epq) + 2*(p+q)
            local_aic(k) = N*log(Epq) + 2*(p + q);

        catch
            local_aic(k) = inf;
        end
    end

    aic_matrix(idx,:) = local_aic;
end

% Select best ARMA(p,q) for each signal
for idx = 1:12
    [~, best_idx] = min(aic_matrix(idx,:));
    best_pq(idx,:) = arma_orders(best_idx,:);
end

% Plot AIC over the 10 candidate models
figure;
t1 = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');

for idx = 1:12
    stim_freq = stim_freq_range(idx);
    aic_values = aic_matrix(idx,:);
    valid_idx = isfinite(aic_values);

    ax = nexttile;
    model_ids = 1:num_models;

    if any(valid_idx)
        plot(ax, model_ids(valid_idx), aic_values(valid_idx), 'b-o', 'LineWidth', 1.2);
        hold(ax, 'on');

        [min_aic, best_idx] = min(aic_values);
        plot(ax, best_idx, min_aic, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');

        p = arma_orders(best_idx,1);
        q = arma_orders(best_idx,2);
        title(ax, ['Freq: ' num2str(stim_freq) ' Hz, ARMA(' num2str(p) ',' num2str(q) ')']);
    else
        title(ax, ['Freq: ' num2str(stim_freq) ' Hz - no valid fit']);
    end

    grid(ax, 'on');
    xlabel(ax, 'Candidate model #');
    ylabel(ax, 'AIC');
    xlim(ax, [1 num_models]);
end

title(t1, 'AIC for 10 Candidate ARMA Models');

%% 7.2) Plot ARMA PSD estimation using the optimal orders

figure;
t2 = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

nfft = 512;

for idx = 1:12
    x = data_filt_all(:,idx);
    stim_freq = stim_freq_range(idx);

    ax = nexttile;

    aic_values = aic_matrix(idx,:);
    valid_idx = isfinite(aic_values);

    if ~any(valid_idx)
        title(ax, [num2str(stim_freq) ' Hz - no valid fit']);
        axis(ax, 'off');
        continue
    end

    [~, best_idx] = min(aic_values);
    p = arma_orders(best_idx,1);    % best p
    q = arma_orders(best_idx,2);    % best q
                 

    try
        Mdl = arima(p,0,q);
        EstMdl = estimate(Mdl, x, 'Display', 'off');

        ar_coeffs = cell2mat(EstMdl.AR);
        ma_coeffs = cell2mat(EstMdl.MA);

        a = [1 ar_coeffs];
        b = [1 ma_coeffs];

        sigma2 = EstMdl.Variance;

        H = freqz(b, a, nfft, Fs);
        Sxx_arma = sigma2 * abs(H).^2;

        f_arma = linspace(0, Fs/2, nfft/2+1);
        Sxx_arma = Sxx_arma(1:nfft/2+1);

        h1 = plot(ax, f_arma, 10*log10(Sxx_arma), 'k-', 'LineWidth', 1.2);
        hold(ax, 'on');
        h2 = xline(ax, stim_freq, 'r--', 'LineWidth', 1);
        h3 = xline(ax, 2*stim_freq, 'b--', 'LineWidth', 1);

        grid(ax, 'on');
        xlim(ax, [0 40]);
        xlabel(ax, 'Freq (Hz)');
        ylabel(ax, 'PSD (dB)');
        title(ax, ['ARMA(' num2str(p) ',' num2str(q) '): ' num2str(stim_freq) ' Hz']);

    catch
        title(ax, [num2str(stim_freq) ' Hz - estimation failed']);
        axis(ax, 'off');
    end
end

title(t2, 'ARMA PSD with Best Model Among 10 Candidates (AIC)');
legend([h1 h2 h3], {'ARMA PSD','1^{st} harmonic','2^{nd} harmonic'}, ...
    'Location','northeast');

%% 8) Welch PSD for all 12 filtered signals

% Welch parameters
winLen = 256;                 % Hamming window length
window = hamming(winLen);
noverlap = round(0.5 * winLen);
nfft = 512;

figure;
t = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

for idx = 1:12
    data = data_sorted(:, idx);
    stim_freq = stim_freq_range(idx);

    % Filter signal
    data_filt = filtfilt(bp, data);

    % Welch PSD
    [pxx, f] = pwelch(data_filt, window, noverlap, nfft, Fs);

    % Plot
    ax = nexttile;
    plot(ax, f, 10*log10(pxx), 'k-', 'LineWidth', 1.2);
    grid(ax, 'on');
    xlim(ax, [0 40]);
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'PSD (dB/Hz)');
    title(ax, ['Welch PSD: ' num2str(stim_freq) ' Hz']);
    xline(ax, stim_freq, 'r--', 'LineWidth', 1);
    xline(ax, 2*stim_freq, 'b--', 'LineWidth', 1);
end

title(t, 'Welch PSD of All 12 Filtered SSVEP Signals');

%% 9) Autocorrelation -> FFT PSD for all 12 filtered signals

%% Autocorrelation -> FFT PSD for all 12 filtered signals

figure;
t = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

for idx = 1:12
    data = data_sorted(:, idx);
    stim_freq = stim_freq_range(idx);

    % Filter signal
    x = filtfilt(bp, data);

    % Remove mean
    x = x - mean(x);

    % Autocorrelation (biased)
    rxx = xcorr(x, 'biased');

    % Use only non-negative lags
    mid = ceil(length(rxx)/2);
    rxx_pos = rxx(mid:end);

    nfft = 500;

    % FFT of autocorrelation
    R = fft(rxx_pos, nfft);

    % PSD estimate
    psd_est = abs(R);

    % Frequency axis
    f = (0:nfft-1)*(Fs/nfft);

    % Keep positive frequencies only
    pos_idx = f >= 0;
    f_pos = f(pos_idx);
    psd_pos = psd_est(pos_idx);

    % Plot
    ax = nexttile;
    plot(ax, f_pos, 10*log10(abs(psd_pos) + eps), 'k-', 'LineWidth', 1.2);
    grid(ax, 'on');
    xlim(ax, [0 40]);
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'PSD (dB)');
    title(ax, ['ACF-FFT PSD: ' num2str(stim_freq) ' Hz']);
    xline(ax, stim_freq, 'r--', 'LineWidth', 1);
    xline(ax, 2*stim_freq, 'b--', 'LineWidth', 1);
end

title(t, 'PSD from Autocorrelation FFT for All 12 Filtered Signals');

%% functions

% Function to compute AIC values for different AR orders
function aic_values = compute_aic_values_AR(data, max_order)
    N = length(data);
    aic_values = zeros(max_order, 1);
    
    for p = 1:max_order
        [a, e] = arburg(data, p);
        aic_values(p) = N * log(e) + 2 * p;
    end
end

