clear 
clc
% close all

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
Is = 2.68*Iwing*L ;         % pitch inertia of store about centre , kg m^2
xf = xs*b ;

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
   
for ii = 1:150
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
    mu = eig(K\[Mt + (0.5*rho*((b^2)/(k^2))*A_hat)]) ;
    % Calculate omega from real part of mu
    omega(:, ii) = sqrt(1./real(mu)) ;
    % Calculate damping from imaginary part of mu
    g(:, ii) = ((omega(:, ii).^2).*imag(mu))/1j ;
    % calculate corresponding speed
    U(:,ii) = (omega(:, ii).*b)/k ;                             % flutter occurs at zero damping
end

% % check for solutions crossing the imaginary axis 
int = find(-0.5*imag(g(2,:))>0) ; 
int = int(end) ;
int = polyfit([U(2,int) U(2,(int+1))], [(-0.5*imag(g(2,int))) (-0.5*imag(g(2,(int+1))))], 1) ;
Uf(int_pos) = abs(int(2))/abs(int(1)) ;

end

% experimental dat afor mass "5"
Span_5 = [0 0.23 0.375 0.54 0.6 0.77 0.85 1] ;
Uf_5 = [101.0412 81.6864 74.0664 71.3232 70.104 73.7616 75.2856 76.80961] ;

% data for mass "7e"
Span_7e = [0 0.23 0.44 0.55 0.7 0.83 1] ;
Uf_7e = [1 0.92 0.87 0.86 0.88 0.94 1.02] ;

% --- Plotting ---
figure
plot(y_L, Uf/Uf(1), 'b-', 'DisplayName', 'Computed (k-method)')
hold on
plot(Span_5, Uf_5/Uf_5(1), '-o', 'DisplayName', 'NACA 1594 (Weight 5)')
xlabel('Dimensionless Spanwise Position (y/L)')
ylabel('Flutter Speed Ratio (U_w / U_0)')
title('Validation against NACA TN 1594')
legend show 
grid on
grid minor

%% Repeat procedure for minimum flutter speed

% Determine the store position for minimum flutter speed
int_pos=find(Uf==min(Uf));
for ii=1:150
    k=ii*0.01;
    Ls=(int_pos-1).*(L/100);
    % Calculate Theodorsen function
    C_theo = besselk(1, (1j*k))./(besselk(0, (1j*k)) + besselk(1, 1j*k)) ;
    % Calculate Deta, Delta, B,C,D and T matrices
    
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
    mu = eig(K\[Mt + (0.5*rho*((b^2)/(k^2))*A_hat)]) ;
    % Calculate omega from real part of mu
    omega(:, ii) = sqrt(1./real(mu)) ;
    % Calculate damping from imaginary part of mu
    g(:, ii) = ((omega(:, ii).^2).*imag(mu))/1j ;
    % calculate corresponding speed
    U(:,ii) = (omega(:, ii).*b)/k ;                             % flutter occurs at zero damping
end

% plots 

figure
hold on 
plot(U(1,:), omega(1,:), '-o')
plot(U(2,:), omega(2,:), '-o')
plot(U(3,:), omega(3,:), '-o')
xlabel('Airspeed (m/s)', 'FontSize', 16)
ylabel('Frequency (Hz)', 'FontSize', 16)
title('Frequency at minumum flutter speed "Mass 5"', 'FontSize', 16, 'FontWeight', 'bold')
legend('mode 1', 'mode 2', 'mode 3')
xlim([0 200])
grid on
grid minor

figure
hold on 
plot(U(1,:), -0.5*imag(g(1,:)), '-o')
plot(U(2,:), -0.5*imag(g(2,:)), '-o')
plot(U(3,:), -0.5*imag(g(3,:)), '-o')
yline(0, '--')
xlabel('Airspeed (m/s)', 'FontSize', 16)
ylabel('-0.5*Damping', 'FontSize', 16)
title('Damping at minumum flutter speed "Mass 5"', 'FontSize', 16, 'FontWeight', 'bold')
legend('mode 1', 'mode 2', 'mode 3')
xlim([0 200])
ylim([-0.1 0.1])
grid on
grid minor