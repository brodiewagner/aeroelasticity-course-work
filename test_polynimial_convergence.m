%% convergence study 
% fix store position at the tip for test
Ls = L; 
xs = 0.625;         
x = xs * b;

% range of polynomial functions 
test_N_values = 2:10; 
num_tests = length(test_N_values);

% arrays to hold the computed frequencies for each test
conv_F1 = zeros(1, num_tests);
conv_F2 = zeros(1, num_tests);
conv_T1 = zeros(1, num_tests);

for k = 1:num_tests
    N = test_N_values(k);
    M = N;                  % using an equal number of bending and torsion modes 
    
    y_L = 0:0.01:1;
    
    % initialize matrices to zero for this iteration size
    Del = zeros(N, N); Dels = zeros(N, N); B = zeros(N, N);
    D = zeros(M, M); D_s = zeros(M, M); T = zeros(M, M);
    C = zeros(N, M); C_s = zeros(N, M);
    
    for i=1:N
        psi_i = (y_L).^(i+1);
        psi_i_store = (Ls/L).^(i+1);
        psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1));
        for j = 1:N
           psi_j = (y_L).^(j+1);
           psi_j_store = (Ls/L).^(j+1);
           psi_j_double = ((j*(j+1))/L.^2).*((y_L).^(j-1));
           
           Del(i,j) = (L)*trapz(y_L, psi_i.*psi_j); 
           Dels(i,j) = psi_i_store.*psi_j_store; 
           B(i,j) = (L)*trapz(y_L, psi_i_double.*psi_j_double); 
        end
    end
    
    for i=1:M
        phi_i = (y_L).^i;
        phi_i_store = (Ls/L).^i;
        phi_i_single = ((i)/(L)).*((y_L).^(i-1));
        for j = 1:M
            phi_j = (y_L).^j;
            phi_j_store = (Ls/L).^j;
            phi_j_single = ((j)/(L)).*((y_L).^(j-1));
            
            D(i,j) = (L)*trapz(y_L, phi_i.*phi_j);
            D_s(i,j) = phi_i_store.*phi_j_store;
            T(i,j) = (L)*trapz(y_L, phi_i_single.*phi_j_single);
        end
    end
    
    for i=1:N
        psi_i = (y_L).^(i+1);
        psi_i_store = (Ls/L).^(i+1);
        for j = 1:M
            phi_j = (y_L).^j;
            phi_j_store = (Ls/L).^j;
            
            C(i,j) = (L)*trapz(y_L, psi_i.*phi_j);
            C_s(i,j) = psi_i_store.*phi_j_store;
        end
    end

    Mwing = [ mbar*Del, -mbar*x_b*C ; -mbar*x_b*(C'), Iwing*D ];
    Mstore = [ ms*Dels, ms*xs*b*C_s ; ms*xs*b*(C_s'), Is*D_s ];
    K = [ EI*B, zeros(N,M) ; zeros(M,N), GJ*T ]; % Note: Fixed bottom-left zero matrix to (M,N)
    
    Mt = Mwing + Mstore;

    [~, Lambda] = eig(K, Mt);
    omega = sqrt(diag(Lambda));
    freqs_hz = sort(omega / (2*pi)); % Sort to grab the lowest 3 modes reliably
    
    conv_F1(k) = freqs_hz(1);
    conv_F2(k) = freqs_hz(2);
    conv_T1(k) = freqs_hz(3);
end

%% Plot Convergence 
figure;
hold on; grid on;

plot(test_N_values, conv_F1, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'w');
plot(test_N_values, conv_F2, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'w');
plot(test_N_values, conv_T1, '-^', 'LineWidth', 2, 'MarkerFaceColor', 'w');

xlabel('Number of Polynomial Functions (N = M)', 'FontSize', 12);
ylabel('Natural Frequency [Hz]', 'FontSize', 12);
title('Convergence Study: Frequency vs. Polynomial Count', 'FontSize', 14, 'FontWeight', 'bold');
legend('Mode 1 (1B)', 'Mode 2 (2B)', 'Mode 3 (1T)', 'Location', 'best');

set(gca, 'FontSize', 12, 'TickDir', 'in');
box on;
hold off;