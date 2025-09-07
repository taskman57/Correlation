function Correlation()
    close all;
    pkg load signal;
    clc;
    % Define RADAR parameters and transmitted signal
    c = 3e8;                            % light speed
    T = 1e-6;                           % Pulse duration
    fs = 50e6;                          % Sampling frequency
    f_c = 1e6;                          % Carrier frequency
    PRI = 50e-6;
    T_period = 1/f_c;
    Np = 1;                            % number of pulses (slow-time)
    Ns = round(PRI*fs);                 % samples per PRI (fast-time)
    N = round(T_period * fs);

    t_pulse = (0:N-1)/fs;
    tx_signal = cos(2*pi*f_c*t_pulse);

    % --- TX matrix 2D: fast-time x slow-time ---
    TX = zeros(Ns, Np);
    for p = 1:Np
        TX(1:N, p) = tx_signal;         % pulse @ start of each PRI
    end

    s1 = subplot(3,1,1);
    plot(TX(:,1));
    set( s1, 'title', 'transmitted signal' , 'fontsize', 20);
    grid minor on;

    % Simulate the echo signal from a target
    target_dist = 1500;                 % Target is 1500m away!

    % 10uS for the pulse to travel to target and back
    time_delay = 2 * target_dist / c;

    % --- RX matrix: delayed copy per pulse (+ noise) ---
    d = round((2*target_dist/c)*fs);    % delay in samples
    RX = zeros(Ns, Np);
    for p = 1:Np
        if d+N-1 <= Ns
            RX(d+1:d+N, p) = tx_signal; % ideal point target, No Doppler yet!
        end
    end

    s2 = subplot(3,1,2);
    plot(RX(:,1));
    grid minor on;
    set(s2, 'title', 'RX (ideal) in PRI #1' , 'fontsize', 20);

    % Add noise to the received signal to make it more realistic
    % noise_std = 0.9;                           % tune this for SNR scenario
    % noisy_rx_sig = 0.05*RX + noise_std*randn(size(RX));

    % --- Add noise based on desired SNR ---
    target_SNR_dB = -40;                                % Examples: -40, -20, -10, 0, +5, +10 dB
    sig_power = mean(abs(RX(:)).^2);                    % Calc average signal power
    target_SNR_lin = 10^(target_SNR_dB/10);             % Conv dB --> linear
    noise_power = sig_power / target_SNR_lin;           % Calc req noise power
    noise_std = sqrt(noise_power);                      % noise standard deviation

    noisy_rx_sig = RX + noise_std*randn(size(RX));      % Apply noise

    s3 = subplot(3,1,3);
    plot(noisy_rx_sig(:,1));
    grid minor on;
    set( s3, 'title', 'RX signal with noise, practical case' , 'fontsize', 20);

    % --- Correlation or Matched filter (per pulse) ---
    % Flip for correlation calculation.
    h = fliplr(tx_signal) / norm(tx_signal);
    MF = zeros(Ns+N-1, Np);
    for p = 1:Np
        MF(:, p) = conv(noisy_rx_sig(:, p), h);         % matched-filter output for pulse p (use noisy RX)
    end

    % --- Range gating & coherent integration (explicit conv-index window) ---
    expected_conv_idx = d + N;                          % 1-based index in conv result
    half_win = round(1.5 * N);                          % small window around mainlobe
    win_idx = max(1, expected_conv_idx - half_win) : min(size(MF,1), expected_conv_idx + half_win);
    coherent_sum = sum(MF(win_idx, :), 2);              % coherent integration across pulses (vector over window)

    % --- Map window indices to range axis for plotting ---
    conv_idx_vec = win_idx;                             % absolute conv indices of the window
    lag_samples_vec = conv_idx_vec - N;                 % lag in samples (0-based)
    lag_axis_s = lag_samples_vec / fs;                  % lag in seconds
    range_axis = (lag_axis_s * c) / 2;                  % range in meters

    % --- plot magnitude and dB with correct x-axis ---
    figure;
    subplot(2,1,1);
    plot(range_axis, abs(coherent_sum), 'LineWidth', 1.2);
    title(sprintf('Coherent sum magnitude (Np=%d)', Np));
    xlabel('Range (m)');
    ylabel('Amplitude');
    grid on;

    subplot(2,1,2);
    plot(range_axis, 20*log10(abs(coherent_sum)+eps), 'LineWidth', 1.2);
    title('Coherent sum (dB)');
    xlabel('Range (m)');
    ylabel('20*log10(|sum|) [dB]');
    grid on;

    % --- compute peak and noise floor (choose regions away from expected peak) ---
    [peak_val, peak_idx] = max(abs(coherent_sum));
    len_cs = length(coherent_sum);
    noise_region = [1:round(0.2*len_cs), round(0.8*len_cs):len_cs];
    noise_rms = sqrt(mean(abs(coherent_sum(noise_region)).^2));

    SNR_linear_amp = peak_val / noise_rms;              % amplitude SNR
    SNR_dB = 20*log10(SNR_linear_amp);
    fprintf('Peak amp = %.3f, noise RMS = %.3f, SNR = %.2f dB\n', peak_val, noise_rms, SNR_dB);

    % --- Range estimate (map conv index back to lag) ---
    peak_conv_idx = win_idx(1) + peak_idx - 1;
    lag_samples_est = peak_conv_idx - N;
    R_est = (lag_samples_est / fs) * c / 2;
    fprintf('Estimated range (coherent): %.2f m\n', R_est);

    % print diagnostics
    fprintf('d = %d, expected_conv_idx = %d, win_idx = [%d %d]\n', d, expected_conv_idx, win_idx(1), win_idx(end));
    fprintf('peak_conv_idx = %d, lag_samples_est = %d\n', peak_conv_idx, lag_samples_est);

end
