% Makale 13
% Signal to Noise Ratio (SNR)
SNR1 = 0:1:30; 

N = 1e6; % # of symbols
A = sqrt(2); % Maximum aplitude of the cos wave being sent

[max_gammak_err01, max_gamma1k_err01, max_gamma2k_err01] = Selection_Scheme_Analytical(SNR1, 2, 0.1);
[max_yk_err01, max_uk_err01] = Selection_Scheme_Simulation(SNR1, A, N, 2, 0.1);

[max_gammak_err10, max_gamma1k_err10, max_gamma2k_err10] = Selection_Scheme_Analytical(SNR1, 2, 10);
[max_yk_err10, max_uk_err10] = Selection_Scheme_Simulation(SNR1, A, N, 2, 10);

semilogy(SNR1, max_gammak_err01, '-* red', SNR1, max_gamma1k_err01, '-x blue', SNR1, max_gamma2k_err01, '-| magenta', SNR1, max_uk_err01, '-o green', SNR1, max_yk_err01, '-diamond cyan');
ylim([10^-4,10^0])
xlim([0,30])
title('BER: BPSK, Comparison of different partial selection schemes at UE = 0.1 and N = 2 in Rayleigh fading channels.')
legend('$$\max{\gamma_k}$$', ...
       '$$\max{\gamma_{1,k}}$$', ...
       '$$\max{\gamma_{2,k}}$$', ...
       '$$\max{\|\mathbf{u}_k\|}$$', ...
       '$$\max{\|\mathbf{y}_k\|}$$', ...
       'Interpreter', 'latex');
xlabel('Es/N0 [dB]');
ylabel('Bit Error Rate (BER)');
grid on

% semilogy(SNR1, max_gammak_err10, '-* red', SNR1, max_gamma1k_err10, '-x blue', SNR1, max_gamma2k_err10, '-| magenta', SNR1, max_uk_err10, '-o green', SNR1, max_yk_err10, '-diamond cyan');
% ylim([10^-4,10^0])
% xlim([0,30])
% title('BER: BPSK, Comparison of different partial selection schemes at UE = 10 and N = 2 in Rayleigh fading channels.')
% legend('$$\max{\gamma_k}$$', ...
%        '$$\max{\gamma_{1,k}}$$', ...
%        '$$\max{\gamma_{2,k}}$$', ...
%        '$$\max{\|\mathbf{u}_k\|}$$', ...
%        '$$\max{\|\mathbf{y}_k\|}$$', ...
%        'Interpreter', 'latex');
% xlabel('Es/N0 [dB]');
% ylabel('Bit Error Rate (BER)');
% grid on

function [max_gammak_err, max_gamma1k_err, max_gamma2k_err] = Selection_Scheme_Analytical(SNR, M_path, UE)
    % Linear SNR
    linear_SNR = 10.^(SNR / 10);
    
    % Error Rate Arrays
    max_gammak_err = zeros(1, length(SNR));
    max_gamma1k_err = zeros(1, length(SNR));
    max_gamma2k_err = zeros(1, length(SNR));
    

    for i = 1:length(SNR)
        % Mean value of SNRs
        gamma1_mean = linear_SNR(i) / (2 * UE + 2);
        gamma2_mean = linear_SNR(i) * UE / (2 * UE + 2);

        sum_gammak = 0;
        sum_gamma1k = 0;
        sum_gamma2k = 0;
        for j = 1:M_path
            sum_gammak = sum_gammak + (nchoosek(M_path, j).* (-1)^j) / sqrt(1 + j/gamma1_mean + j/gamma2_mean);
            sum_gamma1k = sum_gamma1k + (nchoosek(M_path, j).* (-1)^j) / sqrt(1 + j/gamma1_mean + 1/gamma2_mean);
            sum_gamma2k = sum_gamma2k + (nchoosek(M_path, j).* (-1)^j) / sqrt(1 + 1/gamma1_mean + j/gamma2_mean);
        end
        
        max_gammak_err(i) = 1/2 + (1/2) * sum_gammak;
        max_gamma1k_err(i) = 1/2 + (1/2) * sum_gamma1k;
        max_gamma2k_err(i) = 1/2 + (1/2) * sum_gamma2k;
     end
end

function  [max_yk_err, max_uk_err] = Selection_Scheme_Simulation(SNR, A, N, M_path, UE)
    absolute_amplitude = A; % Maximum amplitude of the signal 
    M = 2; % BPSK
    % Energy definitions
    Eb = (absolute_amplitude ^ 2) / 2;
    Es = Eb;

    % Linear SNR
    linear_SNR = 10.^(SNR / 10);

    % Voltage of the sent signals
    voltages = zeros(1,M);
    for i = 1:M
        theta = (i-1) *  180;
      
        voltages(i) = -Es * cosd(theta);
    end

    % Generate random symbols
    s_n = randsrc(1, N, voltages);
    
    % Error Rate Arrays
    max_yk_err = zeros(1, length(SNR));
    max_uk_err = zeros(1, length(SNR));

    % Received Signal Arrays
    x_r = zeros(M_path, N);
    y_d = zeros(M_path, N);
    linked_channel_h = zeros(M_path,N);
    linked_channel_variance = zeros(M_path,N);

    for i = 1:length(SNR)
        % Indirect Paths
        for j = 1:M_path
            h_sr = (1 / sqrt(2)) * (randn(1, N) + 1i * randn(1, N));
            h_rd = (1 / sqrt(2)) * (randn(1, N) + 1i * randn(1, N));
    
            % Power Calculations
            energy_of_received_symbol_sr = mean(abs(h_sr).^2) * Es;
            N0_sr = 2 * energy_of_received_symbol_sr / (linear_SNR(i) / (UE + 1));
            energy_of_received_symbol_rd = mean(abs(h_rd).^2) * Es;
            N0_rd = 2 * energy_of_received_symbol_rd / (linear_SNR(i) * UE / (UE + 1));
            
            % Phase 1 of Indirect Paths
            x_r(j,:) = h_sr .* s_n + sqrt(N0_sr / 2) .* (randn(1,N) + 1i  * randn(1,N));
        
            % Phase 2 of Indirect Paths
            G = sqrt((N0_sr + abs(h_sr).^2));
            s_r = x_r(j,:) ./ G;
            y_d(j,:) = h_rd .* s_r + sqrt(N0_rd / 2) .* (randn(1,N) + 1i  * randn(1,N));

            % Linked Channel Properties
            linked_channel_h(j,:) =  h_sr .* h_rd .* G;
            linked_channel_variance(j,:) = (abs(h_rd).^2 .* N0_rd .* (N0_sr + abs(h_sr).^2).^(-1)) + N0_rd;
        end
        
        % max(|u_k|) Simulation
        [~, u_indices] = max(abs(x_r), [], 1);
        z_n_uk = zeros(1,N);
        decision_u = zeros(1,N);

        for k = 1:N
            z_n_uk(k) =  conj(linked_channel_h(u_indices(k),k)) * y_d(u_indices(k), k);
        end    

        % Decision Circuit
        decision_u(z_n_uk > 0) = voltages(2);
        decision_u(z_n_uk < 0) = voltages(1);

        % SER Simulation
        max_uk_err(i) = sum(decision_u ~= s_n) / N;

        % max(|y_k|) Simulation
        [~, y_indices] = max(abs(y_d), [], 1);
        z_n_yk = zeros(1,N);
        decision_y = zeros(1,N);

        for k = 1:N
             z_n_yk(k) = conj(linked_channel_h(y_indices(k), k)) .* y_d(y_indices(k), k);  
        end

        % Decision Circuit
        decision_y(z_n_yk > 0) = voltages(2);
        decision_y(z_n_yk < 0) = voltages(1);

        % SER Simulation
        max_yk_err(i) = sum(decision_y ~= s_n) / N;

    end     
end