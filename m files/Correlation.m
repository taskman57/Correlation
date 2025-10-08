function Correlation()
    close all;
    pkg load signal;
    clc;

    % =================================================================
    % 1. DEFINE RADAR and TARGET PARAMETERS (PHYSICAL & SYSTEM)
    % =================================================================
    c = 3e8;                                % light speed (m/s)
    T = 1e-6;                               % Pulse duration (s)
    fs = 50e6;                              % Sampling frequency (Hz)
    f_c = 10e9;                             % Carrier frequency (10 GHz - X-band)

    % --- MEDIUM RANGE PARAMETERS ---
    PRI = 50e-6;                            % PRI for max range 7.5 km
    target_dist = 7000;                     % Target at 7030 m (Multiple of 3m res.)

    % --- System & Target Physical Parameters (RADAR Range Equation) ---
    % Pt adjusted for a Single-Pulse SNR_in ~ -17 dB at 7 km
    Pt = 100.00;                            % FINAL Pt: 160 W, Realistic Medium Power
    G_dB = 20;                              % Antenna Gain (dBi)
    G = 10^(G_dB / 10);                     % Antenna Gain (unitless)
    lambda = c / f_c;                       % Wavelength (m)
    sigma = 1;                              % Target RCS (m^2)
    B = 1/T;                                % Receiver Bandwidth (Hz, Matched Filter)
    F = 2;                                  % Noise Figure (Unitless, approx 3 dB)
    L = 5;                                  % System Loss Factor (Unitless, approx 7 dB)

    % Constants
    k = 1.380649e-23;                       % Boltzmann's Constant (J/K)
    T0 = 290;                               % Standard Noise Temperature (K)

    % Integration
    Np = 8;                                 % <<< Set to 1 for UNDETECTABLE, then 100 for SUCCESS >>>

    % --- Practical Noise Components ---
    f_ripple = 120;
    V_ripple_rms = 5e-7;
    R_ref = 50;
    P_ripple = V_ripple_rms^2 / R_ref;

    % =================================================================
    % 2. TX SIGNAL GENERATION (COMPLEX BASEBAND)
    % =================================================================
    Ns = round(PRI*fs);                     % Samples per PRI (2500 samples)
    N = round(T * fs);                      % Pulse samples (N=50, MF Gain is 16.99 dB)
    t_pulse = (0:N-1)/fs;
    B_chirp = 20e6;                         % Chirp bandwidth (20 MHz ??? ~7.5 m range resolution)
    k_chirp = B_chirp / T;                  % Chirp rate
    tx_signal_bb = exp(1i * pi * k_chirp * 1 * t_pulse.^2);   % LFM chirp

    TX = zeros(Ns, Np);
    for p = 1:Np
        TX(1:N, p) = tx_signal_bb;
    end

    s1 = subplot(3,1,1);
    plot(real(TX(:,1)));
    set( s1, 'title', 'Transmitted Signal (Real Part of Baseband)' , 'fontsize', 14);
    grid minor on;

    % =================================================================
    % 3. RX SIGNAL GENERATION & SCALING
    % =================================================================

    time_delay = 2 * target_dist / c;
    d = round(time_delay * fs);

    % --- Calculate Single-Pulse SNR using RADAR Range Equation ---
    R = target_dist;
    Denominator_Loss = (4 * pi)^3 * R^4 * L;
    S_watts = (Pt * G^2 * lambda^2 * sigma) / Denominator_Loss;

    N_thermal_watts = k * T0 * B * F * L;

    SNR_single_ratio = S_watts / N_thermal_watts;
    SNR_single_dB = 10 * log10(SNR_single_ratio);

    fprintf(' Gain of TX/RX antenna and Calculated System Parameters (R=%.0fm):\n', R);
    fprintf(' Antenna gain for both TX & RX: %d dBi\n', G_dB);
    fprintf('  - Pt: %.2f W, MF Gain: %.2f dB\n', Pt, 10*log10(N));
    fprintf('  - Single-Pulse SNR (Theoretical): %.2f dB\n', SNR_single_dB);
    fprintf('  - S_watts: %.2e W, N_thermal: %.2e W\n', S_watts, N_thermal_watts);

    % --- Calculate Noise Standard Deviations ---
    sigma_n_thermal = sqrt(N_thermal_watts / 2);
    sigma_n_ripple = sqrt(P_ripple / 2);

    % --- RX matrix: delayed copy per pulse (Complex) ---
    RX_ideal = zeros(Ns, Np); % <<< FIXED: Initialized RX_ideal here >>>
    for p = 1:Np
        if d+N-1 <= Ns
            RX_ideal(d+1:d+N, p) = tx_signal_bb;
        end
    end

    % --- Scale the IDEAL RX signal to match the required S_watts ---
    % Note: RX_ideal is now defined from line 84
    P_ideal_current = mean(abs(RX_ideal(:)).^2);
    scale_factor = sqrt(S_watts / P_ideal_current);
    RX_scaled = RX_ideal * scale_factor;

    s2 = subplot(3,1,2);
    plot(real(RX_scaled(:,1)));
    grid minor on;
    set(s2, 'title', sprintf('RX (Scaled to S=%.2eW) in PRI #1 (Visually undetectable)', S_watts) , 'fontsize', 14);


    % 4. Add AWGN and Ripple Noise

    noise_thermal = sigma_n_thermal * (randn(size(RX_scaled)) + 1i * randn(size(RX_scaled)));

    t_full = (0:Ns*Np-1)/fs;
    ripple_i = sigma_n_ripple * sin(2*pi*f_ripple*t_full);
    ripple_q = sigma_n_ripple * cos(2*pi*f_ripple*t_full);
    noise_ripple_full = (ripple_i + 1i * ripple_q).';
    noise_ripple = reshape(noise_ripple_full, Ns, Np);

    noisy_rx_sig = RX_scaled + noise_thermal + noise_ripple;

    s3 = subplot(3,1,3);
    plot(real(noisy_rx_sig(:,1)));
    grid minor on;
    set( s3, 'title', sprintf('Noisy RX (SNR=%.2fdB)', SNR_single_dB) , 'fontsize', 14);

    % =================================================================
    % 5. MATCHED FILTERING & COHERENT INTEGRATION
    % =================================================================

    h = conj(fliplr(tx_signal_bb)) / norm(tx_signal_bb);
    MF = zeros(Ns+N-1, Np);
    for p = 1:Np
        MF(:, p) = conv(noisy_rx_sig(:, p), h);
    end

    expected_conv_idx = d + N;
    half_win = round(2.5 * N);
    win_idx = max(1, expected_conv_idx - half_win) : min(size(MF,1), expected_conv_idx + half_win);
    coherent_sum = sum(MF(win_idx, :), 2);

    lag_samples_vec = win_idx - N;
    lag_axis_s = lag_samples_vec / fs;
    range_axis = (lag_axis_s * c) / 2;

    single_pulse_result = MF(win_idx, 1);

    % =================================================================
    % 6. PLOTTING & DIAGNOSTICS
    % =================================================================

    figure('Name', sprintf('Integration Comparison (Np=%d, SNR=%.2fdB)', Np, SNR_single_dB));

    s_comp = subplot(2,1,1);
    plot(range_axis, abs(single_pulse_result), 'b', 'LineWidth', 1.2);
    hold on;
    plot(range_axis, abs(coherent_sum), 'r', 'LineWidth', 1.2);
    title(sprintf('Coherent Integration vs. Single Pulse (Np=%d)', Np));
    xlabel('Range (m)');
    ylabel('Amplitude');
    legend('Single Pulse', sprintf('Coherent Sum (Gain: %.1fdB)', 10*log10(Np)));
    grid on;

    s_db = subplot(2,1,2);
    plot(range_axis, 20*log10(abs(single_pulse_result)+eps), 'b', 'LineWidth', 1.2);
    hold on;
    plot(range_axis, 20*log10(abs(coherent_sum)+eps), 'r', 'LineWidth', 1.2);
    title('Matched Filter Output (dB)');
    xlabel('Range (m)');
    ylabel('20*log10(|R|) [dB]');
    legend('Single Pulse', sprintf('Coherent Sum (Np=%d)', Np));
    grid on;

    % --- compute observed SNR ---
    [peak_val_single, ~] = max(abs(single_pulse_result));
    [peak_val_coherent, peak_idx] = max(abs(coherent_sum));

    len_cs = length(coherent_sum);
    noise_region = 1 : max(1, peak_idx - round(3*N));

    if isempty(noise_region) || length(noise_region) < 10
        noise_region = [1:round(0.2*len_cs), round(0.8*len_cs):len_cs];
    end

    noise_rms_coherent = sqrt(mean(abs(coherent_sum(noise_region)).^2));

    SNR_dB_observed = 20*log10(peak_val_coherent / noise_rms_coherent);

    noise_rms_single = noise_rms_coherent / sqrt(Np);
    SNR_dB_single_observed = 20*log10(peak_val_single / noise_rms_single);

    fprintf('\n--- Observed Results ---\n');
    fprintf('MF Output SNR (Theoretical, Np=1): %.2f dB\n', SNR_single_dB + 10*log10(N));
    fprintf('Coherent Integration Gain (Theoretical): +%.2f dB\n', 10*log10(Np));
    fprintf('Single Pulse Observed SNR (MF Output): %.2f dB\n', SNR_dB_single_observed);
    fprintf('Coherent Sum Observed SNR (MF Output): %.2f dB\n', SNR_dB_observed);
    fprintf('Observed SNR Improvement (Coherent): %.2f dB\n', SNR_dB_observed - SNR_dB_single_observed);

    % --- Range estimate ---
    peak_conv_idx = win_idx(1) + peak_idx - 1;
    lag_samples_est = peak_conv_idx - N;
    R_est = (lag_samples_est / fs) * c / 2;
    fprintf('Estimated range: %.2f m (Target: %.0f m)\n', R_est, target_dist);

end
