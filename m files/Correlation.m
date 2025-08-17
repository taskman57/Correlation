function Correlation()
    close all;
    pkg load signal;
    clc;
    % Define RADAR parameters and transmitted signal
    c = 3e8;                            % light speed
    T = 1e-6;                           % Pulse duration
    fs = 10e6;                          % Sampling frequency
    t_pulse = 0 : 1/fs : T;
    tx_signal = cos(2*pi*5e6*t_pulse);  % TX signal (5 MHz carrier)
    s1 = subplot(3,1,1);
    plot(tx_signal);
    set( s1, 'title', 'transmitted signal' , 'fontsize', 20);
    grid minor on;
    % Simulate the echo signal from a target
    target_dist = 1500;                 % Target is 1500m away!

    % 10uS for the pulse to travel to target and back
    time_delay = 2 * target_dist / c;

    % Received signal, transmitted and back
    t_sim = 0 : 1/fs : 2*time_delay;    % Longer time vector for simulation
    rx_signal = zeros(size(t_sim));
    % Find the start index for the delayed signal
    delay_samples = round(time_delay * fs);
    % Place the delayed pulse into the received signal array
    rx_signal(delay_samples + 1 : delay_samples + length(tx_signal)) = tx_signal;
    s2 = subplot(3,1,2);
    plot(rx_signal);
    grid minor on;
    set( s2, 'title', 'RX from TX signal, no noise!' , 'fontsize', 20);

    % Add noise to the received signal to make it more realistic
    noise = 0.9 * randn(size(rx_signal));
    noisy_rx_sig = rx_signal + 0.5*noise;
    s3 = subplot(3,1,3);
    plot(noisy_rx_sig);
    grid minor on;
    set( s3, 'title', 'RX signal with noise, practical case' , 'fontsize', 20);

    % Perform Cross-Correlation to find the target's echo
    % xcorr is used for simulation.
    [correlation_result, lags] = xcorr(noisy_rx_sig, tx_signal);

    % Find the peak of the correlation result
    % The peak tells us the time delay
    [max_corr_value, peak_index] = max(abs(correlation_result));
    peak_lag = lags(peak_index);

    % target range calculation
    % Convert the lag (in samples) back to time delay, and then to distance.
    measured_time_delay = peak_lag / fs;
    measured_range = (measured_time_delay * c) / 2;

    % Print the results
    fprintf(' Simulated Target Range: %.2f meters\n', target_dist);
    fprintf(' Calculated Range via Correlation: %.2f meters\n', measured_range);

    % Plotting the results to visualize
    figure;
    plot(lags/fs, correlation_result);
    title('Cross-Correlation Result', 'fontsize', 20);
    xlabel('Time Lag (s)', 'fontsize', 20);
    ylabel('Correlation', 'fontsize', 20);
    grid minor on;
    % Uncomment above lines to visualize
end
