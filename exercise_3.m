clear 
clc
close all

%% Define Parameters
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

mbar = mw/L ;               % mass per unit span of wing kg/m
I_total = 4.349e-3;         % Total wing pitch inertia
Iwing = I_total / L;        % Moment of inertia per unit span
rho = mbar/(pi*b*b*32.6);   % air density (kg/m^3)
e = (EA/c)-(1/4);           % elastic axis local relative quarter chord
a = (e*c)/b-0.5;            % elastic axis position

% parameters for wing store 5 
ms = 0.636*mw ;             % store mass, kg
xs = -0.687 ;               % position of store ahead of Elastic Axis, m
Is = 2.68*Iwing*L ;         % pitch inertia of store about centre , kg·m^2
xf = xs*b ;
weight_no = '5' ;

% % parameters for wing store 7e 
% ms =  0.954*mw ;          % store mass, kg
% xs = 0.034 ;              % position of store ahead of Elastic Axis, m
% Is = 1.56*Iwing*L ;       % pitch inertia of store about centre, kg·m^2
% xf = xs*b ;
% weight_no = '7e' ;

%% Assumed Shapes Method
% store position at 100 points along span
for int_pos = 1:1:101
    Ls = (int_pos-1).*(L/100) ;
    y_L = 0:0.01:1 ;
    N = 3 ;
    M = 2 ;
    
    % calculate Delta, Deltas, B, C, D and T matrices setup functions for bending, 3 bending and 2 torsion modes 
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

Mstore = [    ms*Dels,     ms*xs*b*C_s   ; 
           ms*xs*b*(C_s'),   Is*D_s   ]   ;

%Determince expression for the stiffness matrix, unaffected by the addition of a store
K = [    EI*B,      zeros(N,M)  ;
     (zeros(N,M))',    GJ*T   ] ;

% Total mass matrix achieved by combining the mass matrices of store and wing
    Mt = Mwing + Mstore ;
   
% Divergence and Flutter (k-method)
for ii = 1:1:150
    k = ii*0.01;
    % Calculate Theodorsen function
    C_theo = besselk(1,(1j*k))./(besselk(0,(1j*k)) + besselk(1,1j*k)) ;
    % Calculate the A_mat, B_Mat, C_mat matrices
    A_mat = 2*pi*b*(k^2)*[    Del           a*b*C             ;
                            a*b*(C') (b^2)*((a^2) + 1/8)*D  ] ;
    B_mat = -2*pi*k*1j*[       2*C_theo*Del             -b*(1+2*(0.5-a)*C_theo)*C       ; 
                         2*b*(0.5+a)*C_theo*(C') (b^2)*(0.5-a)*(1-2*(0.5+a)*C_theo)*D ] ;
    C_mat = -2*pi*b*[ zeros(N,N)      -2*C_theo*C        ; 
                      zeros(N,M)' -b*(1+(2*a))*C_theo*D ] ;
    A_hat = A_mat + B_mat + C_mat ;
    mu = eig(K\(Mt + (0.5*rho*((b^2)/(k^2))*A_hat))) ;
    % Calculate omega from real part of mu
    omega(:, ii) = sqrt(1./real(mu)) ;
    % Calculate damping from imaginary part of mu
    g(:, ii) = ((omega(:, ii).^2).*imag(mu))/1j ;
    % calculate corresponding speed
    U(:,ii) = (omega(:, ii).*b)/k ;                             % flutter occurs at zero damping
end

% check for solutions crossing the imaginary axis 
int = find(-0.5*imag(g(2,:))>0) ; 
int = int(end) ;
int = polyfit([U(2,int) U(2,(int+1))], [(-0.5*imag(g(2,int))) (-0.5*imag(g(2,(int+1))))], 1) ;
Uf(int_pos) = abs(int(2))/abs(int(1)) ;

end

%% Plotting & Experimental Data from NACA TN 1594

% experimental dat afor mass "5"
    Span_5 = [0 0.234043 0.382979 0.553191 0.617021 0.787234 0.87234 1] ;
    Uf_5_fps = [331.5 288 243 234.5 230 242 247 252] ;
    % Uf_5 = [101.0412 81.6864 74.0664 71.3232 70.104 73.7616 75.2856 76.80961] ;        % From graph 
% convert feet per second to meters per second
    fps_to_mps = 3.28083989 ;
    Uf_5 = Uf_5_fps*fps_to_mps ;

% data for mass "7e"
    % From graph
        % Span_7e = [0, 0.23093922, 0.43976074, 0.54911888, 0.7296852, 0.83078590, 1] ;
        % Uf_7e = [1, 0.934953, 0.854259, 0.8897015, 0.886371, 0.9493911, 1.023524] ;
    % From table data
        Span_7e = [0 0.234043 0.446809 0.617021 0.744681 0.851064 1.021277] ;
        Uf_7e_fps = [331.5 308 283 294 292 310 337] ;
        Uf_7e = Uf_7e_fps*fps_to_mps ;

% Plot
figure('Name', sprintf('Validated Flutter Ratio - Weight %s', weight_no))
plot(y_L, Uf/Uf(1), '-o', 'DisplayName', 'Computed (k-method)', 'LineWidth', 2)         
hold on
plot(Span_5, Uf_5/Uf_5(1), '^', ...                                             % plotting experimental data
     'DisplayName', 'Exp. Data - Weight 5', ...
     'MarkerEdgeColor', [0.87 0.33 0], ...
     'LineWidth', 1.5, 'MarkerSize', 10)
p = polyfit(Span_5, Uf_5/Uf_5(1), 4); Uf_5_bestfit = polyval(p, y_L);           % set up for best fit curve
plot(y_L, Uf_5_bestfit, '-', ...                                                % plotting best fit curve
    'LineWidth', 1.75, ...
    'DisplayName', 'Best Fit Curve - Weight 5', ...
    'Color', [0.87 0.33 0])
xlabel('y/L', 'FontSize', 15)
ylabel('$\biggl(\frac{U_w}{U_{F0}}\biggr)$', ...
       'Interpreter', 'latex', 'FontSize', 15, 'Rotation', 0)
title('Validation against NACA TN 1594 - Weight 5')
set(gca, 'FontName', 'helvetica')
legend('FontSize', 10)
grid on
grid minor

% figure('Name', sprintf('Validated Flutter Ratio - Weight %s', weight_no))
% plot(y_L, Uf/Uf(1), '-', 'DisplayName', 'Computed (k-method)', 'LineWidth', 2)
% hold on
% plot(Span_7e, Uf_7e/Uf_7e(1), '^', ...
%      'DisplayName', 'Exp. Data - Weight 7', ...
%      'MarkerEdgeColor', [0.87 0.33 0], ...
%      'LineWidth', 1.5, 'MarkerSize', 10)
% p = polyfit(Span_7e, Uf_7e/Uf_7e(1), 5); Uf_7e_bestfit = polyval(p, y_L);
% plot(y_L, Uf_7e_bestfit, '-', ...
%     'LineWidth', 1.75, ...
%     'DisplayName', 'Best Fit Curve - Weight 5', ...
%     'Color', [0.87 0.33 0])
% xlabel('y/L', 'FontSize', 15)
% ylabel('$\biggl(\frac{U_w}{U_{F0}}\biggr)$', ...
%        'Interpreter', 'latex', 'FontSize', 15, 'Rotation', 0)
% title('Validation against NACA TN 1594 - Weight 7e')
% set(gca, 'FontName', 'helvetica')
% xlim([0 1])
% legend('FontSize', 10)
% grid on
% grid minor

%% Repeat procedure for minimum flutter speed

% Determine the store position for minimum flutter speed
int_pos=find(Uf==min(Uf));
for ii=1:1:150
    k=ii*0.01;
    Ls=(int_pos-1).*(L/100);
    % Calculate Theodorsen function
    C_theo = besselk(1, (1j*k))./(besselk(0, (1j*k)) + besselk(1, 1j*k)) ;

    % set up functions for bending 
    y_L = 0:0.01:1 ;
    N = 3 ;
    M = 2 ;
    
    % calculate Delta, Deltas, B, C, D and T matrices setup functions for bending, 3 bending and 2 torsion modes 
    for i=1:N
        psi_i = (y_L).^(i+1) ;                                      % ith bending function for wing
        psi_i_store = (Ls/L).^(i+1) ;                               % ith bending funtion for store
        psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second derivative (psi'')
        for j = 1:N
           psi_j = (y_L).^(j+1) ;                                   % jth bedning function of wing
           psi_j_store = (Ls/L).^(j+1) ;                            % jth bending mode of store
           psi_j_double = ((j*(j+1))/L.^2).*((y_L).^(j-1)) ;        % jth wing bending function second derivative (psi'') 
           % Del, Dels and B matrices
           Del(i,j) = (L-0)*trapz(y_L,psi_i.*psi_j) ;                   % Delta for wing
           Dels(i,j) = psi_i_store.*psi_j_store ;                       % Delta for store
           B(i,j) = (L-0)*trapz(y_L,psi_i_double.*psi_j_double) ;   
        end
    end
    
    % setup i and j torsion for wing and store, 2 torsion modes
    for i=1:M
        phi_i = (y_L).^i ;                                          % ith torsional function for wing
        phi_i_store = (Ls/L).^i ;                                   % ith torsional function for store
        phi_i_single = ((i)/(L)).*((y_L).^(i-1)) ;                  % ith torsional function first derivative (phi')
        for j = 1:M
            phi_j = (y_L).^j ;                                      % jth torsional function for wing 
            phi_j_store = (Ls/L).^j ;                               % jth torsional function for store
            phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;              % jth torsional function first derivative (phi')
            % D, D_store and T matrices
            D(i,j) = (L-0)*trapz(y_L, phi_i.*phi_j) ;                   % D matrix for wing
            D_s(i,j) = phi_i_store.*phi_j_store ;                       % D matrix for store
            T(i,j) = (L-0)*trapz(y_L, phi_i_single.*phi_j_single) ;     % T matrix
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

Mstore = [    ms*Dels,     ms*xs*b*C_s   ; 
           ms*xs*b*(C_s'),   Is*D_s   ]   ;

% Determine expression for the stiffness matrix, unaffected by the addition of a store
K = [    EI*B,      zeros(N,M)  ;
     (zeros(N,M))',    GJ*T   ] ;

% Total mass matrix achieved by combining the mass matrices of store and wing
    Mt = Mwing + Mstore ;          

    % Calculate the A_mat, B_Mat, C_mat matrices
    A_mat = 2*pi*b*(k^2)*[    Del           a*b*C             ;
                            a*b*(C') (b^2)*((a^2) + 1/8)*D  ] ;
    B_mat = -2*pi*k*1j*[       2*C_theo*Del             -b*(1+2*(0.5-a)*C_theo)*C       ; 
                         2*b*(0.5+a)*C_theo*(C') (b^2)*(0.5-a)*(1-2*(0.5+a)*C_theo)*D ] ;
    C_mat = -2*pi*b*[ zeros(N,N)      -2*C_theo*C        ; 
                      zeros(N,M)' -b*(1+(2*a))*C_theo*D ] ;
    A_hat = A_mat + B_mat + C_mat ;
    mu = eig(K\(Mt + (0.5*rho*((b^2)/(k^2))*A_hat))) ;
    % Calculate omega from real part of mu
    omega(:, ii) = sqrt(1./real(mu)) ;
    % Calculate damping from imaginary part of mu
    g(:, ii) = ((omega(:, ii).^2).*imag(mu))/1j ;
    % calculate corresponding speed
    U(:,ii) = (omega(:, ii).*b)/k ;                             % flutter occurs at zero damping
end


% plots 
figure('Name', sprintf('Minimum Flutter - Weight %s', weight_no))
tiledlayout(1, 2)
nexttile
% figure('Name', 'Frequency at Min. Flutter')
hold on 
plot(U(1,:), omega(1,:), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 1')
plot(U(2,:), omega(2,:), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 2')
plot(U(3,:), omega(3,:), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 3')
xlabel('Airspeed (m/s)', 'FontSize', 15)
ylabel('Frequency (Hz)', 'FontSize', 15)
title('Frequency at minumum flutter speed' , 'FontSize', 15, 'FontWeight', 'bold')
legend('Location', 'southeast')
xlim([0 200])
grid on
grid minor
nexttile
% figure('Name', 'Damping at Min. Flutter')
hold on 
plot(U(1,:), -0.5*imag(g(1,:)), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 1')
plot(U(2,:), -0.5*imag(g(2,:)), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 2')
plot(U(3,:), -0.5*imag(g(3,:)), '-o', 'LineWidth', 1.5, 'DisplayName', 'Mode 3')
yline(0, '-.', 'HandleVisibility', 'off')
xlabel('Airspeed (m/s)', 'FontSize', 15)
ylabel('-0.5*Damping', 'FontSize', 15)
title('Damping at minumum flutter speed', 'FontSize', 15, 'FontWeight', 'bold')
legend('Location', 'southeast')
xlim([0 200])
ylim([-0.1 0.1])
grid on
grid minor