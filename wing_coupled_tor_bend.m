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

% store position at 100 points along span
 for int_pos = 1:1:101
    Ls = (int_pos-1).*(L/100) ;
    % calculate Delta, Deltas, B, C, D and T matrices.
    % setup functions for bending, 3 bending and 2 torsion modes
    y_L = 0:0.1:1 ;
    N = 5 ;                                                      % assume as a first test the bending 
    M = 6 ;                                                      % and torsion modes

    for i=1:N
        psi_i = (y_L).^(i+1) ;                                      % ith bending function for wing
        psi_i_store = (Ls/L).^(i+1) ;                               % ith bending funtion for store
        psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second derivative (psi'')
        for j = 1:N
           psi_j = (y_L).^(j+1) ;                                   % jth bending function of wing
           psi_j_store = (Ls/L).^(j+1) ;                            % jth bending mode of store
           psi_j_double = ((j*(j+1))/L.^2).*((y_L).^(j-1)) ;        % jth wing bending function second derivative (psi'') 
           % Del, Dels and B matrices
           Del(i,j) = (L-0)*trapz(y_L,psi_i.*psi_j) ;              % Delta for wing
           Dels(i,j) = psi_i_store.*psi_j_store ;                  % Delta for store
           B(i,j) = (L-0)*trapz(y_L,psi_i_double.*psi_j_double) ;  % 
        end
    end
    % setup i and j torsion for wing and store, 2 torsion modes
    for i=1:M
        phi_i = (y_L).^i ;                                          % ith torsional function for wing
        phi_i_store = (Ls/L).^i ;                                   % ith torsional function for store
        phi_i_single = ((i)/(L)).*((y_L).^(i-1)) ;                  % ith torsional function first derivative (phi')
        for j = 1:M
            phi_j = (y_L).^j ;                                          % jth torsional function for wing 
            phi_j_store = (Ls/L).^j ;                                   % jth torsional function for store
            phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;                  % jth torsional function first derivative (phi')
            % D, D_store and T matrices
            D(i,j) = (L-0)*trapz(y_L, phi_i.*phi_j) ;                       % D matrix for wing
            D_s(i,j) = phi_i_store.*phi_j_store ;                           % D matrix for store
            T(i,j) = (L-0)*trapz(y_L, phi_i_single.*phi_j_single) ;         % T matrix
        end
    end
    % Define functions for i bending and j torsion to determine C matrix
    for i=1:N
        psi_i = (y_L).^(i+1) ;                                      % ith bending function for wing
        psi_i_store = (Ls/L).^(i+1) ;                               % ith bending funtion for store
        psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second derivative (psi'')
        for j = 1:M
            phi_j = (y_L).^j ;                                      % jth torsional function for wing
            phi_j_store = (Ls/L).^j ;                               % jth torsional function for store 
            phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;              % jth torsional function first derivative (phi')
            % C and C_store matrices
            C(i,j) = (L-0)*trapz(y_L, psi_i.*phi_j) ;
            C_s(i,j) = psi_i_store.*phi_j_store ;
        end
    end

% Expressions for wing and store mass matrix
Mwing = [   mbar*Del,   -mbar*x_b*C   ;
         -mbar*x_b*(C'),  Iwing*D   ] ;

Mstore = [    ms*Dels,     ms*x*C_s   ; 
           ms*x*(C_s'),   Is*D_s   ]   ;

%Determince expression for the stiffness matrix, unaffected by the addition of a store
K = [    EI*B,      zeros(N,M)  ;
     (zeros(N,M))',    GJ*T   ] ;

% Total mass matrix achieved by combining the mass matrices of store and wing
    Mt = Mwing + Mstore ;

% Determine natural frequencies
    [V,Lambda] = eig(K,Mt) ;
    
    % We only want to save the eigenvectors at 11 specific points along the span.
    if int_pos == 1
        Vout(:,:,1) = V ;                   % Save the Ls = 0 (clean wing) case as index 1
        
    elseif mod(int_pos-1, 10) == 0
        store_idx = (int_pos-1)/10 + 1 ;    % For int_pos 11, 21, 31... calculate the corresponding index (2 through 11)
        Vout(:,:,store_idx) = V ;           % Save the 5x5 eigenvector matrix into the 3D Vout array
    end

    omega = sqrt(diag(Lambda)) ;            % rad/s frequency
    f(:,int_pos) = omega/(2*pi) ;           % transforms the Frequency from radians to Hz
end
    % Data from NACA TN1594 for plotting
  
    %Deflection and twist of the elastic axis create 11 store positions along span
    y_L=0:0.1:1;
    % select mode
    mode = 2 ;
    % select store position
    store_pos = 2 ;
    EA_z = zeros(1,length(y_L)) ;
    theta_EA = zeros(1,length(y_L)) ;
    % Determine the amount of bending
    for j = 1:N
        a = 0.01*Vout(j,mode,store_pos) ;
        psi_j = (y_L).^(j+1) ;
        EA_z = EA_z + psi_j*a ;
    end
    % Determine the amount of twist
    for j = N+1:(M+N)
        b = 0.01*Vout(j,mode,store_pos) ;
        phi_j = (y_L).^(j-N) ;
        theta_EA = theta_EA+(phi_j*b)*(180/pi) ;
    end

% Geometry --------------------------------------------------------------------------------------------------------------

    % Determine leading edge and trailing edge x and z coordinates
    Lead_x = 0.1 ;
    Lead_z= 0 ;
    Trail_x = -0.1 ; 
    Trail_z = 0 ;
    Store_x = 0 ;
    Store_z = 0 ;
       
    % Calculate physical span locations
    y_span = L * y_L; 

    % Calculate the deformed coordinates for the whole wing (Using cosd and sind because theta_EA is in degrees)
    x_LE = Lead_x * cosd(theta_EA);
    z_LE = Lead_x * sind(theta_EA) + EA_z;
    
    x_TE = Trail_x * cosd(theta_EA);
    z_TE = Trail_x * sind(theta_EA) + EA_z;
    
    x_EA = zeros(1, length(y_L));
    z_EA = EA_z;

%% Experimental Data
%          y/l      f_b1   f_b2   f_t
data = [
    0,     0,        6.45, 39.2   47.3    ;
    0,     0,        6.43, 39.2   39.2    ;
    11,    0.234043, 6.64, 30.41, 40      ;
    11,    0.234043, 6.64, 30.15, 40.21   ;
    20.5,  0.43617,  6.14, 35.6,  35.6    ; 
    20.5,  0.43617,  6.15, 35.6,  35.6    ;
    26,    0.553191, 5.56, 33.69, 24.27   ; 
    26,    0.553191, 5.54, 34,    22.91   ;
    29,    0.617021, 5.22, 35.1,  24.74   ;
    29,    0.617021, 5.3,  36.11, 24.63   ;
    41,    0.87234,  4.17, 38.5,  23.9    ;
    41,    0.87234,  4.17, 38.3,  23.9    ;
    44,    0.93617,  3.86, 35.83, 23.59   ;
    44,    0.93617,  3.86, 36,    23.59   ;
    44,    0.93617,  3.83, 35.8,  23.38   ;
    44,    0.93617,  3.8,  25.64, 23.33   ;
    44,    0.93617,  3.84, 35.83, 23.29   ;
    44,    0.93617,  3.87, 36.35, 23.52   ;
    47,    1,        3.65, 34.35, 22.5    ;
    47,    1,        3.68, 34.25, 22.5    ;
    47,    1,        3.59, 33.13, 20.63   ;
    47,    1,        3.62, 33.74, 22.6    ;
    47,    1,        3.59, 33.74, 22.02   ;
    47,    1,        3.61, 34.09, 22.83 ] ;

store_pos_in = data(:, 1);
store_pos_yL = data(:, 2);
exp_B1 = data(:, 3);            % experimental data for first bending 
exp_B2 = data(:, 4);            % experimental data for second bending 
exp_T  = data(:, 5);            % experimental data for torsion

% Find the unique y/L positions
[unique_yL, ~, idx] = unique(store_pos_yL);

% Calculate the Mean for each unique position
mean_B1 = accumarray(idx, exp_B1, [], @mean);
mean_B2 = accumarray(idx, exp_B2, [], @mean);
mean_T  = accumarray(idx, exp_T,  [], @mean);

% Calculate the Standard Deviation (Error) for each unique position
std_B1 = accumarray(idx, exp_B1, [], @std);
std_B2 = accumarray(idx, exp_B2, [], @std);
std_T  = accumarray(idx, exp_T,  [], @std);

%% 3D PLOT --------------------------------------------------------------------------------------------------------------

    figure('name', '3D visualisation of Wing');
    hold on;
    grid on;
    view(-35, 300) ;                  % Sets the 3D viewing angle to match your reference image
    
    % rigid wing 
    plot3([Lead_x Lead_x], [0 L], [0 0], '--k');                    % Leading Edge
    plot3([Trail_x Trail_x], [0 L], [0 0], '--k');                  % Trailing Edge
    plot3([0 0], [0 L], [0 0], '--k');                              % Elastic Axis
    
    % plot chordwise ribs
    for i = 1:length(y_L)
        plot3([Lead_x, Trail_x], [y_span(i), y_span(i)], [0, 0], '--k');
    end

    % deformed wing 
    plot3(x_LE, y_span, z_LE, '-k'); % Deformed Leading Edge
    plot3(x_TE, y_span, z_TE, '-k'); % Deformed Trailing Edge
    plot3(x_EA, y_span, z_EA, '-k'); % Deformed Elastic Axis
    
    % plot chordwise ribs
    for i = 1:length(y_L)
        plot3([x_LE(i), x_TE(i)], [y_span(i), y_span(i)], [z_LE(i), z_TE(i)], '-k');
    end

    % text denoting leading and trailing edges
    text(x_LE(6), y_span(6), z_LE(6) + 0.005, 'LE', 'FontSize', 14, 'FontWeight', 'bold');    % Index 6 corresponds to 50% 
    text(x_TE(6), y_span(6), z_TE(6) + 0.005, 'TE', 'FontSize', 14, 'FontWeight', 'bold');    % span (y/L = 0.5)

    % store marker
    y_store = y_span(store_pos);                                    % Calculate store 
    x_store_rot = x * cosd(theta_EA(store_pos));                    % location on the 
    z_store = x * sind(theta_EA(store_pos)) + EA_z(store_pos);      % deformed wing
    
    % plot the store as a large square 's' marker
    plot3(x_store_rot, y_store, z_store, 'sk', 'MarkerSize', 15, 'LineWidth', 2);

    % figure formatting 
    
    % force desired view angle 
    view(130, 25) ;
    
    % axis labels
    xlabel('x (m)', 'FontSize', 14);
    ylabel('Span (m)', 'FontSize', 14);
    zlabel('z (m)', 'FontSize', 14);

    % axis limits 
    xlim([-0.2, 0.2]);
    ylim([0, 1.4]);
    
    % make the 3D bounding box rectangular
    pbaspect([1 2.5 0.8]);              % [x-width, y-length, z-height]. This makes Span 2.5x visually longer than Chord.
    
    % title based on selected mode and store position
    title_str = sprintf('Mode %d for store position (y/L)=%.1f', mode, y_L(store_pos));
    title(title_str, 'FontSize', 16, 'FontWeight', 'bold');
   
    hold off;

%% 2D Frequency plot 
% 
% % sort frequencies at every spanwise station from lowest to highest
% f_sorted = sort(f, 1);          % ensures F1 is always the lowest line, F2 is the middle, etc
% 
% % create x-axis array for computed lines 
% x_span = 0:0.01:1;
% 
% figure('Name', 'Nat. Freq. for first 3 modes');
% hold on;
% 
% % plot the computed lines (rows 1, 2, and 3 of f_sorted)
% plot(x_span, f_sorted(1,:), '-', 'LineWidth', 2.5, 'DisplayName', 'Computed B1');  
% plot(x_span, f_sorted(2,:), '-', 'LineWidth', 2.5, 'DisplayName', 'Computed B2');  
% plot(x_span, f_sorted(3,:), '-', 'LineWidth', 2.5, 'DisplayName', 'Computed T1');  
% 
% % Plot the experimental data markers with ERROR BARS
% % Syntax: errorbar(x, y, error, 'LineStyle', 'none', 'Marker', ...)
% errorbar(unique_yL, mean_B1, std_B1, '^k', 'MarkerSize', 15, 'LineWidth', 1, 'DisplayName', 'Experimental B1');
% errorbar(unique_yL, mean_B2, std_B2, 'ok', 'MarkerSize', 15, 'LineWidth', 1, 'DisplayName', 'Experimental B2');
% errorbar(unique_yL, mean_T, std_T, 'sk', 'MarkerSize', 15, 'LineWidth', 1, 'DisplayName', 'Experimental T1');
% 
% % figure formatting 
% xlabel('Store position (y/L)', 'FontSize', 16);
% ylabel('Frequency [Hz]', 'FontSize', 16);
% 
% % title based on selected N and M values
% title_str = sprintf('First 3 modal frequencies vs Store Position (N = %d, M = %d)', N, M) ;
% title(title_str, 'FontSize', 18, 'FontWeight', 'bold');
% 
% % axis limits
% xlim([-0.01 1]);
% ylim([0 50]);
% 
% grid on
% grid minor
% legend('Location', 'northeast', 'FontSize', 15, "Position", [0.8388 0.7589 0.1555 0.2167]);
% 
% % clean up the axes 
% box off;
% set(gca, 'FontSize', 20, 'TickDir', 'in');
% hold off;

