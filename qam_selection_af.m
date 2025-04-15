% Zhao, Y., Adve, R., & Lim, T. J. (2006). 
% Symbol error rate of selection amplify-and-forward relay systems. 
% IEEE Transactions on Communications, 10(11), 757-759.
% Fig.1 & Fig.2

% Signal to Noise Ratio (SNR)
SNR1 = 5:1:25; 

N = 1e6; % # of symbols
A = sqrt(2); % Maximum aplitude of the cos wave being sent

ser_approx1 = Approximated_Results(SNR1, 4, 2);
ser_approx2 = Approximated_Results(SNR1, 16, 2);
ser_analytical1 = Theoretical_Results(SNR1, A, N, 4, 2);
ser_analytical2 = Theoretical_Results(SNR1, A, N, 16, 2);

semilogy(SNR1, ser_approx1, '-x blue', SNR1, ser_approx2, '-x blue', SNR1, ser_analytical1, '-o red', SNR1, ser_analytical2, '-o red');
ylim([10^-8,1.2])
xlim([5,25])
title('SER: QAM, SER of S-AF Network with Two Relays')
legend('Approximation (M = 4)', 'Approximation (M = 16)', 'Theoretical (M = 4)', 'Theoretical (M = 16)');
xlabel('Es/N0 [dB]');
ylabel('Bit Error Rate (BER)');
grid on

function [ser_analytical] = Theoretical_Results(SNR, A, N, M, M_path)

    % Constants regarding the M-QAM Modulations 
    if M == 4
        c = 2;
    elseif M == 16
        c = 0.8;
    end
    
    linear_SNR = 10.^(SNR / 10) / M_path; % linearized SNR

    ser_analytical = zeros(1, length(SNR));

    Q_func = @(x) 0.5 * erfc(x / sqrt(2));

    PDF_gamma_r = @(gamma_r, gamma_val) ...
    (M_path) * (gamma_r.^(M_path)) ./ (gamma_val.^(M_path + 1)) .* exp(-gamma_r ./ gamma_val);

    for i = 1:length(SNR)

        instant_snr = linear_SNR(i);
        integrand = @(gamma_r) Q_func(sqrt(c * gamma_r)) .* PDF_gamma_r(gamma_r, instant_snr);
            
        ser_analytical(i) = integral(integrand, 0, Inf);
    end
       
end

function [ser_approx] = Approximated_Results(SNR, M, M_path)

    linear_SNR = 10.^(SNR / 10) / M_path; % Linear SNR
    
    % Constants regarding the M-QAM Modulations 
    % Q((2i-1) * sqrt(3*Eb*log_2(M)/(N0*(M-1))))
    if M == 4
        c = 2;
    elseif M == 16
        c = 0.8;
    end

    ser_approx = zeros(1, length(SNR));
    lambda_i = sqrt(0.5);
    ksi_i = sqrt(0.5);
    lambda0 = sqrt(0.5);
    for i = 1:length(SNR)
        product = 1;
        for j = 1: M_path
            product = product *(lambda_i+ksi_i);
        end
        Pe = factorial(2*M_path + 1) .* lambda0 .* product / (factorial(M_path+1) * (2*c*linear_SNR(i))^(M_path + 1));
        
        ser_approx(i) = Pe;
    end


end
