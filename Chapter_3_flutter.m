clear 
clc
% close all

%% initialise parameters
% wing parameters
mw = 3.48*0.4535924 ;       % wing mass (kg)
AR = 6 ;                    % Aspect Ratio
c = 8*0.0254 ;              % root chord (m)
b = c/2 ;                   % root semi-chord (m)
L = c*AR ;                  % wing semi-span (m)
EA = 0.437*c ;              % non dimensional elastic axis position chords
IA = 0.454*c ;              % non dimensional inertial axis position chords
x_b = IA - EA ; 

lb_in_to_Nm2 = 4.44822*(0.0254^2) ;     % 1 lb = 4.44822 N, 1 in = 0.0254 m
EI = 140700*lb_in_to_Nm2 ;              % bending stiffness (N·m^2)
GJ = 69200*lb_in_to_Nm2 ;               % torsional stiffness (N·m^2/rad)
mbar = mw/L ;                           % mass per unit span of wing kg/m
I_total = 4.349e-3;                     % total wing pitch inertia
Iwing = I_total / L;                    % moment of inertia per unit span
rho = mbar/(pi*b*b*32.6);               % air density (kg/m^3)
e = (EA/c)-(1/4);                       % elastic axis local relative quarter chord
a = (e*c)/b-0.5;                        % elastic axis position

%% initialise experimental data

fps_to_mps = 3.28083989 ;

% Weight 5 Data
Span_5 = [0 0.234043 0.382979 0.553191 0.617021 0.787234 0.87234 1] ;
Uf_5_fps = [331.5 288 243 234.5 230 242 247 252] ;
Uf_5 = Uf_5_fps * fps_to_mps ;

% Weight 7e Data 
Span_7e = [0 0.234043 0.446809 0.617021 0.744681 0.851064 1.021277] ;
Uf_7e_fps = [331.5 308 283 294 292 310 337] ;
Uf_7e = Uf_7e_fps * fps_to_mps ;

% discretise span y/L 
y_L = 0:0.01:1 ;

%% Loop through both configurations
for config = 1:2
    
    if config == 1
        % Parameters for wing store 5 
        ms = 0.636*mw ;             % store mass, kg
        xs = -0.687 ;               % position of store ahead of Elastic Axis, m
        Is = 2.68*Iwing*L ;         % pitch inertia of store about centre , kg·m^2
        xf = xs*b ;
        weight_no = '5' ;
    else
        % Parameters for wing store 7e 
        ms =  0.954*mw ;          % store mass, kg
        xs = 0.034 ;              % position of store ahead of Elastic Axis, m
        Is = 1.56*Iwing*L ;       % pitch inertia of store about centre, kg·m^2
        xf = xs*b ;
        weight_no = '7e' ;
    end
        
    %% assumed shapes method
    % store position at 100 points along span
    for int_pos = 1:1:101
        Ls = (int_pos-1).*(L/100) ;
        N = 5 ;
        M = 6 ;

        % calculate Delta, Deltas, B, C, D and T matrices setup functions for bending, N bending and M torsion modes 
        for i=1:N
            psi_i = (y_L).^(i+1) ;                                      % ith bending function for wing
            psi_i_store = (Ls/L).^(i+1) ;                               % ith bending funtion for store
            psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second derivative (psi'')
            for j = 1:N
               psi_j = (y_L).^(j+1) ;                                   % jth bedning function of wing
               psi_j_store = (Ls/L).^(j+1) ;                            % jth bending mode of store
               psi_j_double = ((j*(j+1))/L.^2).*((y_L).^(j-1)) ;        % jth wing bending function second derivative (psi'') 
               % Del, Dels and B matrices
               Del(i,j) = (L-0)*trapz(y_L,psi_i.*psi_j) ;               % Delta for wing
               Dels(i,j) = psi_i_store.*psi_j_store ;                   % Delta for store
               B(i,j) = (L-0)*trapz(y_L,psi_i_double.*psi_j_double) ;   % B matrix
            end
        end
        
        % setup i and j torsion for wing and store, 2 torsion modes
        for i=1:M
            phi_i = (y_L).^i ;                                          % ith torsional function for wing
            phi_i_store = (Ls/L).^i ;                                   % ith torsional function for store
            phi_i_single = ((i)/(L)).*((y_L).^(i-1)) ;                   % ith torsional function first derivative (phi')
            for j = 1:M
                phi_j = (y_L).^j ;                                      % jth torsional function for wing 
                phi_j_store = (Ls/L).^j ;                               % jth torsional function for store
                phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;               % jth torsional function first derivative (phi')
                % D, D_store and T matrices
                D(i,j) = (L-0)*trapz(y_L, phi_i.*phi_j) ;                % D matrix for wing
                D_s(i,j) = phi_i_store.*phi_j_store ;                    % D matrix for store
                T(i,j) = (L-0)*trapz(y_L, phi_i_single.*phi_j_single) ;       % T matrix
            end
        end
        
        % define functions for i bending and j torsion to determine C matrix
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
        
        % expressions for wing and store mass matrix
        Mwing = [   mbar*Del,   -mbar*x_b*C   ;
                 -mbar*x_b*(C'),  Iwing*D   ] ;
        Mstore = [    ms*Dels,     ms*xs*b*C_s   ; 
                   ms*xs*b*(C_s'),   Is*D_s   ]   ;

        % determine expression for the stiffness matrix, unaffected by the addition of a store
        K = [    EI*B,      zeros(N,M)  ;
             (zeros(N,M))',    GJ*T   ] ;
        
        % total mass matrix achieved by combining the mass matrices of store and wing
        Mt = Mwing + Mstore ;
       
        % divergence and flutter (k-method)
        for ii = 1:1:150
            k = ii*0.01;
            % calculate Theodorsen function
            C_theo = besselk(1,(1j*k))./(besselk(0,(1j*k)) + besselk(1,1j*k)) ;
            % calculate the A_mat, B_Mat, C_mat matrices
            A_mat = 2*pi*b*(k^2)*[    Del           a*b*C             ;
                                    a*b*(C') (b^2)*((a^2) + 1/8)*D  ] ;
            B_mat = -2*pi*k*1j*[       2*C_theo*Del             -b*(1+2*(0.5-a)*C_theo)*C       ; 
                                 2*b*(0.5+a)*C_theo*(C') (b^2)*(0.5-a)*(1-2*(0.5+a)*C_theo)*D ] ;
            C_mat = -2*pi*b*[ zeros(N,N)      -2*C_theo*C        ; 
                              zeros(N,M)' -b*(1+(2*a))*C_theo*D ] ;
            A_hat = A_mat + B_mat + C_mat ;
            mu = eig(K\(Mt + (0.5*rho*((b^2)/(k^2))*A_hat))) ;
            
            % calculate omega from real part of mu
            omega(:, ii) = sqrt(1./real(mu)) ;
            
            % calculate damping from imaginary part of mu
            g(:, ii) = ((omega(:, ii).^2).*imag(mu))/1j ;
            
            % calculate corresponding speed
            U_speed(:,ii) = (omega(:, ii).*b)/k ;                             % flutter occurs at zero damping
        end
        
        % check for solutions crossing the imaginary axis 
        int_cross = find(-0.5*imag(g(2,:))>0) ; 
        int_cross = int_cross(end) ;
        interp_fit = polyfit([U_speed(2,int_cross) U_speed(2,(int_cross+1))], [(-0.5*imag(g(2,int_cross))) (-0.5*imag(g(2,(int_cross+1))))], 1) ;
        Uf(int_pos) = abs(interp_fit(2))/abs(interp_fit(1)) ;
    end
    
    % Store the flutter curve for the combined plot later
    Uf_results(config, :) = Uf;
    
    %% repeat for minimum flutter speed plots
    % determine the store position for minimum flutter speed
    int_pos=find(Uf==min(Uf));
    for ii=1:1:150
        k=ii*0.01;
        Ls=(int_pos-1).*(L/100);
        C_theo = besselk(1, (1j*k))./(besselk(0, (1j*k)) + besselk(1, 1j*k)) ;
        
        for i=1:N
            psi_i = (y_L).^(i+1) ;                                      
            psi_i_store = (Ls/L).^(i+1) ;                               
            psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           
            for j = 1:N
               psi_j = (y_L).^(j+1) ;                                   
               psi_j_store = (Ls/L).^(j+1) ;                            
               psi_j_double = ((j*(j+1))/L.^2).*((y_L).^(j-1)) ;        
               Del(i,j) = (L-0)*trapz(y_L,psi_i.*psi_j) ;                   
               Dels(i,j) = psi_i_store.*psi_j_store ;                       
               B(i,j) = (L-0)*trapz(y_L,psi_i_double.*psi_j_double) ;   
            end
        end
        
        for i=1:M
            phi_i = (y_L).^i ;                                          
            phi_i_store = (Ls/L).^i ;                                   
            phi_i_single = ((i)/(L)).*((y_L).^(i-1)) ;                  
            for j = 1:M
                phi_j = (y_L).^j ;                                      
                phi_j_store = (Ls/L).^j ;                               
                phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;              
                D(i,j) = (L-0)*trapz(y_L, phi_i.*phi_j) ;                   
                D_s(i,j) = phi_i_store.*phi_j_store ;                       
                T(i,j) = (L-0)*trapz(y_L, phi_i_single.*phi_j_single) ;     
            end
        end
        
        for i=1:N
            psi_i = (y_L).^(i+1) ;                                      
            psi_i_store = (Ls/L).^(i+1) ;                               
            psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           
            for j = 1:M
                phi_j = (y_L).^j ;                                      
                phi_j_store = (Ls/L).^j ;                               
                phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;              
                C(i,j) = (L-0)*trapz(y_L, psi_i.*phi_j) ;
                C_s(i,j) = psi_i_store.*phi_j_store ;
            end
        end

        Mwing = [ mbar*Del, -mbar*x_b*C; 
                 -mbar*x_b*(C'), Iwing*D ] ;
        Mstore = [ ms*Dels, ms*xs*b*C_s; 
                   ms*xs*b*(C_s'), Is*D_s ] ;
        K = [ EI*B, zeros(N,M); 
            (zeros(N,M))', GJ*T ] ;
        Mt = Mwing + Mstore ;          
        A_mat = 2*pi*b*(k^2)*[ Del, a*b*C; a*b*(C'), (b^2)*((a^2) + 1/8)*D ] ;
        B_mat = -2*pi*k*1j*[ 2*C_theo*Del, -b*(1+2*(0.5-a)*C_theo)*C; 
                             2*b*(0.5+a)*C_theo*(C'), (b^2)*(0.5-a)*(1-2*(0.5+a)*C_theo)*D ] ;
        C_mat = -2*pi*b*[ zeros(N,N), -2*C_theo*C; 
                          zeros(N,M)', -b*(1+(2*a))*C_theo*D ] ;
        A_hat = A_mat + B_mat + C_mat ;
        mu = eig(K\(Mt + (0.5*rho*((b^2)/(k^2))*A_hat))) ;
        omega(:, ii) = sqrt(1./real(mu)) ;
        g(:, ii) = ((omega(:, ii).^2).*imag(mu))/1j ;
        U_speed(:,ii) = (omega(:, ii).*b)/k ;                             
    end
    
    % calculate the exact y/L ratio for the minimum flutter speed
    min_y_L = (int_pos - 1) / 100; 

    % plots for individual minimum flutter
    figure('Name', sprintf('Minimum Flutter - Weight %s', weight_no))
    t = tiledlayout(1, 2); % Assign tiledlayout to a variable
    
    % add main title across entire figure
    title(t, sprintf('Critical Flutter Position for Weight %s: y/L = %.2f', weight_no, min_y_L), ...
          'FontSize', 16, 'FontWeight', 'bold')
    nexttile
    hold on 
    plot(U_speed(1,:), omega(1,:), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 1')
    plot(U_speed(2,:), omega(2,:), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 2')
    plot(U_speed(3,:), omega(3,:), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 3')
    xlabel('Airspeed, U_\infty (m/s)', 'FontSize', 20)
    ylabel('Frequency, f (Hz)', 'FontSize', 20)
    title(sprintf('Frequency at min flutter (Store %s)', weight_no), 'FontSize', 15, 'FontWeight', 'bold')
    legend('Location', 'southeast', 'FontSize', 15)
    xlim([0 200])
    grid on
    grid minor
    nexttile
    hold on 
    plot(U_speed(1,:), -0.5*imag(g(1,:)), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 1')
    plot(U_speed(2,:), -0.5*imag(g(2,:)), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 2')
    plot(U_speed(3,:), -0.5*imag(g(3,:)), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 3')
    yline(0, '-.', 'HandleVisibility', 'off')
    xlabel('Airspeed, U_\infty (m/s)', 'FontSize', 20)
    ylabel('-0.5*Damping', 'FontSize', 20)
    title(sprintf('Damping at min. flutter (Store %s)', weight_no), 'FontSize', 15, 'FontWeight', 'bold')
    legend('Location', 'southeast', 'FontSize', 15)
    xlim([0 150])
    ylim([-0.1 0.1])
    grid on
    grid minor
    
end % End of loop over configurations

%% Combined Plotting: Validated Flutter Ratio - Weights 5 and 7e
figure('Name', 'Validated Flutter Ratio - Stores 5 and 7e', 'Position', [100, 100, 800, 600])

% Define colors to distinguish the weights
color_5 = [0 0.427 0.831] ;       % MATLAB Blue
color_7e = [0.286 0.678 0];       % MATLAB Green

% Plot Computed K-Method Results
Uf_5_comp = Uf_results(1,:);
Uf_7e_comp = Uf_results(2,:);

plot(y_L, Uf_5_comp/Uf_5_comp(1), '-', ...
    'Color', color_5, ...
    'DisplayName', 'Numerical (Store 5)', ...
    'LineWidth', 2)         
hold on
plot(y_L, Uf_7e_comp/Uf_7e_comp(1), '-', ...
    'Color', color_7e, ...
    'DisplayName', 'Numerical (Store 7e)', ...
    'LineWidth', 2)

% Plot Experimental Data
plot(Span_5, Uf_5/Uf_5(1), '^', ...
     'DisplayName', 'Exp. Data (Store 5)', ...
     'MarkerEdgeColor', color_5, ...
     'MarkerFaceColor', 'none', ...
     'LineWidth', 1.5, 'MarkerSize', 8)
 
plot(Span_7e, Uf_7e/Uf_7e(1), 'square', ...
     'DisplayName', 'Exp. Data (Store 7e)', ...
     'MarkerEdgeColor', color_7e, ...
     'MarkerFaceColor', 'none', ...
     'LineWidth', 1.5, 'MarkerSize', 8)

% Plot Best Fit Curves
p5 = polyfit(Span_5, Uf_5/Uf_5(1), 5); 
Uf_5_bestfit = polyval(p5, y_L);           
plot(y_L, Uf_5_bestfit, '--', ...
    'Color', color_5, ...
    'LineWidth', 1.5, ...
    'HandleVisibility', 'off') 

p7e = polyfit(Span_7e, Uf_7e/Uf_7e(1), 5); 
Uf_7e_bestfit = polyval(p7e, y_L);
plot(y_L, Uf_7e_bestfit, '--', ...
    'Color', color_7e, ...
    'LineWidth', 1.5, ...
    'HandleVisibility', 'off')

% Formatting
xlabel('y/L', 'FontSize', 15)
ylabel('$\biggl(\frac{U_f}{U_{fr}}\biggr)$', ...
       'Interpreter', 'latex', ...
       'FontSize', 18, ...
       'Rotation', 0, ...
       'HorizontalAlignment', 'right')
title('Validation against NACA TN 1594 - Stores 5 & 7e', 'FontSize', 14)
set(gca, 'FontName', 'helvetica')
legend('FontSize', 11, 'Location', 'best')
xlim([0 1])
grid on
grid minor
hold off