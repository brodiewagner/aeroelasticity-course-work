clc
clear
close all

%%
% wing parameters
mw = 1.5818 ;                   % mass of the wing
c = 0.2032 ;                    % chord length 
b = c/2 ;                       % half chord
L = 1.2192 ;                    % length of the wing
EA = 0.437*c ;                  % position of elastic axis
IA = 0.454*c ;                  % position of inertial axis
x_b = IA-EA ;                   % 
EI = 404.76 ;                   % flexural rigidity of wing (Nm^2)
GJ = 199.076 ;                  % torsional rigidity of wing (Nm^2/rad)
mbar = mw/L ;                   % mass per unit length of wing (kg/m)
Iwing = 4.349e-3/L ;            % wing pitch moment of inertia (kg m^2/m)
rho = mbar/(pi*b*b*32.6) ;      

% parameters for wing store 4 
ms = 0.636*mw ;                    % store 4 mass, kg
xs = 0.625 ;                       % position of store ahead of Elastic Axis, m
Is = 1.91*Iwing*L ;                % pitch inertia of store about centre , kg m^2
x = xs*b ;

%% senstivity analysis
disp('Running 2D Sensitivity Analysis (Varying N and M)...')

% 1. Set our test condition: Store at 11 inches (11/47 * L = 0.234043 * L)
Ls_test = 0.234043 * L; 
y_L_int = 0:0.01:1; % Integration vector for trapz

% Experimental Data from NACA TN 1594 (Table II, Run 9, Weight 4 at 11 inches)
f1B_exp = 6.41;  % 1st Bending [Hz]
f2B_exp = 30.95; % 2nd Bending [Hz]
f1T_exp = 40.24; % 1st Torsion [Hz]

% 2. We will test N and M from 1 to 8
N_values = 1:8;
M_values = 1:8;

% Matrices to store results (rows = N values, cols = M values)
freq_results_1 = zeros(length(N_values), length(M_values));
freq_results_2 = zeros(length(N_values), length(M_values));
freq_results_3 = zeros(length(N_values), length(M_values));

% Nested loops to evaluate every combination of N and M
for n_idx = 1:length(N_values)
    for m_idx = 1:length(M_values)
        N_test = N_values(n_idx);
        M_test = M_values(m_idx);
        
        % Initialize matrices for this loop
        Del = zeros(N_test,N_test); Dels = zeros(N_test,N_test); B = zeros(N_test,N_test);
        D = zeros(M_test,M_test); D_s = zeros(M_test,M_test); T = zeros(M_test,M_test);
        C = zeros(N_test,M_test); C_s = zeros(N_test,M_test);
        
        % Build Bending Matrices
        for i=1:N_test
            psi_i = (y_L_int).^(i+1); 
            psi_i_store = (Ls_test/L).^(i+1); 
            psi_i_double = ((i*(i+1))/L^2).*((y_L_int).^(i-1)); 
            for j=1:N_test
                psi_j = (y_L_int).^(j+1); 
                psi_j_store = (Ls_test/L).^(j+1); 
                psi_j_double = ((j*(j+1))/L^2).*((y_L_int).^(j-1)); 
                
                Del(i,j) = L*trapz(y_L_int, psi_i.*psi_j); 
                Dels(i,j) = psi_i_store.*psi_j_store; 
                B(i,j) = L*trapz(y_L_int, psi_i_double.*psi_j_double); 
            end
        end
        
        % Build Torsion Matrices
        for i=1:M_test
            phi_i = (y_L_int).^i; 
            phi_i_store = (Ls_test/L).^i; 
            phi_i_single = (i/L).*((y_L_int).^(i-1)); 
            for j=1:M_test
                phi_j = (y_L_int).^j; 
                phi_j_store = (Ls_test/L).^j; 
                phi_j_single = (j/L).*((y_L_int).^(j-1)); 
                
                D(i,j) = L*trapz(y_L_int, phi_i.*phi_j); 
                D_s(i,j) = phi_i_store.*phi_j_store; 
                T(i,j) = L*trapz(y_L_int, phi_i_single.*phi_j_single); 
            end
        end
        
        % Build Coupled Matrices
        for i=1:N_test
            psi_i = (y_L_int).^(i+1); 
            psi_i_store = (Ls_test/L).^(i+1); 
            for j=1:M_test
                phi_j = (y_L_int).^j; 
                phi_j_store = (Ls_test/L).^j; 
                
                C(i,j) = L*trapz(y_L_int, psi_i.*phi_j);
                C_s(i,j) = psi_i_store.*phi_j_store;
            end
        end
        
        % Assemble Global Matrices
        Mwing_test = [mbar*Del, -mbar*x_b*C; -mbar*x_b*(C'), Iwing*D];
        Mstore_test = [ms*Dels, ms*x*C_s; ms*x*(C_s'), Is*D_s];
        Mt_test = Mwing_test + Mstore_test;
        K_test = [EI*B, zeros(N_test,M_test); zeros(M_test,N_test), GJ*T];
        
        % Solve Eigenvalue Problem
        [~, Lambda_test] = eig(K_test, Mt_test);
        f_hz = real(sqrt(diag(Lambda_test))) / (2*pi); % Keep real part for stability
        f_sorted = sort(f_hz); 
        
        % Safety check: Pad f_sorted with NaNs if it has fewer than 3 elements
        % (This happens when N+M < 3, e.g., N=1 and M=1)
        f_padded = [f_sorted; NaN(max(0, 3 - length(f_sorted)), 1)];
        
        % Store first 3 modes safely
        freq_results_1(n_idx, m_idx) = f_padded(1);
        freq_results_2(n_idx, m_idx) = f_padded(2);
        freq_results_3(n_idx, m_idx) = f_padded(3);
    end
end

%% 3. Plot the 2D Sensitivity Results
[M_grid, N_grid] = meshgrid(M_values, N_values);

% % --- Subplot 1: First Mode ---
% figure('Name', 'N and M Convergence vs Experimental Data', 'Position', [100, 100, 1400, 450]);
% surf(M_grid, N_grid, freq_results_1, 'FaceAlpha', 0.8); 
% hold on;
% % Plot Experimental Plane
% surf(M_grid, N_grid, f1B_exp*ones(size(freq_results_1)), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
% title('Computed Mode 1 vs. Exp 1st Bending', 'FontSize', 12);
% xlabel('Torsional Polynomials (M)'); ylabel('Bending Polynomials (N)'); zlabel('Frequency [Hz]');
% view(-45, 30); grid on;

% --- Subplot 2: Second Mode ---
figure('Name', 'N and M Convergence vs Experimental Data', 'Position', [100, 100, 1400, 450]);
surf(M_grid, N_grid, freq_results_2, 'FaceAlpha', 0.75); 
hold on;
% Plot Experimental Plane
surf(M_grid, N_grid, f2B_exp*ones(size(freq_results_2)), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('Computed 2nd Bending Mode Frequency Convergence', 'FontSize', 20);
xlabel('Torsional Polynomials (M)', 'FontSize', 15, 'Rotation', 30); 
ylabel('Bending Polynomials (N)', 'FontSize', 15, 'Rotation', -30); 
zlabel('Frequency [Hz]', 'FontSize', 15);
view([133.5 25.0]);
legend('Convergence', 'Experimental Data', 'fontsize', 15)
set(gca, 'FontSize', 20);
grid on;

% % --- Subplot 3: Third Mode ---
figure('Name', 'N and M Convergence vs Experimental Data', 'Position', [100, 100, 1400, 450]);
surf(M_grid, N_grid, freq_results_3, 'FaceAlpha', 0.75); 
hold on;
% Plot Experimental Plane
surf(M_grid, N_grid, f1T_exp*ones(size(freq_results_3)), 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('Computed 1st Torsional Mode Frequency Convergence', 'FontSize', 20);
xlabel('Torsional Polynomials (M)', 'FontSize', 25, 'Rotation', 30); 
ylabel('Bending Polynomials (N)', 'FontSize', 25, 'Rotation', -30); 
zlabel('Frequency [Hz]', 'FontSize', 25);
view([133.5 25.0]);
legend('Convergence', 'Experimental Data', 'fontsize', 15)
set(gca, 'FontSize', 20);
grid on;

% sgtitle('Rayleigh-Ritz Convergence Sensitivity (Store Mass 4 @ 0.234 L)', 'FontSize', 16, 'FontWeight', 'bold');