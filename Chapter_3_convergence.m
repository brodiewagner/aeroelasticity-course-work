clear 
clc
close all

%% define Parameters
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

%% compute structural matrices for store at y/L = 0.61
Ls = 0.61 * L ; 
y_L = 0:0.01:1 ;
N = 5 ; 
M = 6 ; 

for i=1:N
    psi_i = (y_L).^(i+1) ;
    psi_i_store = (Ls/L).^(i+1) ;
    psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;
    for j = 1:N
       psi_j = (y_L).^(j+1) ;
       psi_j_store = (Ls/L).^(j+1) ;
       psi_j_double = ((j*(j+1))/L.^2).*((y_L).^(j-1)) ;
       Del(i,j) = L*trapz(y_L, psi_i.*psi_j) ;
       Dels(i,j) = psi_i_store.*psi_j_store ;
       B(i,j) = L*trapz(y_L, psi_i_double.*psi_j_double) ;
    end
end

for i=1:M
    phi_i = (y_L).^i ;
    phi_i_store = (Ls/L).^i ;
    phi_i_single = (i/L).*((y_L).^(i-1)) ;
    for j = 1:M
        phi_j = (y_L).^j ;
        phi_j_store = (Ls/L).^j ;
        phi_j_single = (j/L).*((y_L).^(j-1)) ;
        D(i,j) = L*trapz(y_L, phi_i.*phi_j) ;
        D_s(i,j) = phi_i_store.*phi_j_store ;
        T(i,j) = L*trapz(y_L, phi_i_single.*phi_j_single) ;
    end
end

for i=1:N
    psi_i = (y_L).^(i+1) ;
    psi_i_store = (Ls/L).^(i+1) ;
    for j = 1:M
        phi_j = (y_L).^j ;
        phi_j_store = (Ls/L).^j ;
        C(i,j) = L*trapz(y_L, psi_i.*phi_j) ;
        C_s(i,j) = psi_i_store.*phi_j_store ;
    end
end

Mwing = [   mbar*Del,   -mbar*x_b*C   ;
         -mbar*x_b*(C'),  Iwing*D   ] ;
Mstore = [    ms*Dels,     ms*xs*b*C_s   ; 
           ms*xs*b*(C_s'),   Is*D_s   ]   ;
Mt = Mwing + Mstore ;
K = [    EI*B,      zeros(N,M)  ;
     zeros(M,N),    GJ*T   ] ;

%% convergence Study for k-resolution

k_max = 1.5; % Define maximum reduced frequency

% array of different resolutions (number of k points)
k_resolutions = [10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 1500]; 
Uf_converged = zeros(size(k_resolutions));
step_sizes = zeros(size(k_resolutions));

for i = 1:length(k_resolutions)
    num_points = k_resolutions(i);
    dk = k_max / num_points; % calculate the step size
    step_sizes(i) = dk;
    
    for ii = 1:num_points
        k = ii * dk;
        
        C_theo = besselk(1,(1j*k)) / (besselk(0,(1j*k)) + besselk(1,1j*k)) ;
        
        A_mat = 2*pi*b*(k^2)*[    Del           a*b*C             ;
                                a*b*(C') (b^2)*((a^2) + 1/8)*D  ] ;
        B_mat = -2*pi*k*1j*[       2*C_theo*Del             -b*(1+2*(0.5-a)*C_theo)*C       ; 
                             2*b*(0.5+a)*C_theo*(C') (b^2)*(0.5-a)*(1-2*(0.5+a)*C_theo)*D ] ;
        C_mat = -2*pi*b*[ zeros(N,N)      -2*C_theo*C        ; 
                          zeros(N,M)' -b*(1+(2*a))*C_theo*D ] ;
                          
        A_hat = A_mat + B_mat + C_mat ;
        mu = eig(K\(Mt + (0.5*rho*((b^2)/(k^2))*A_hat))) ;
        
        % Store Results
        omega(:, ii) = sqrt(1./real(mu)) ;
        g(:, ii) = ((omega(:, ii).^2).*imag(mu))/1j ;
        U(:,ii) = (omega(:, ii).*b)/k ; 
    end
    
    % find where damping crosses 0 (flutter condition) for Mode 2
    damp_mode2 = -0.5 * imag(g(2,:));
    int_idx = find(damp_mode2 > 0);
    
    if ~isempty(int_idx) && int_idx(end) < num_points
        int_idx = int_idx(end); % Last index before crossing
        
        % linear interpolation for a more precise zero crossing
        U_val1 = U(2, int_idx);
        U_val2 = U(2, int_idx+1);
        g_val1 = damp_mode2(int_idx);
        g_val2 = damp_mode2(int_idx+1);
        
        poly_fit = polyfit([U_val1 U_val2], [g_val1 g_val2], 1);
        Uf_converged(i) = abs(poly_fit(2)) / abs(poly_fit(1));
    else
        Uf_converged(i) = NaN; % Flutter did not occur or out of bounds
    end
end

%% plots
figure('Name', 'k-Resolution Convergence Study');
tiledlayout(1,2, 'TileSpacing', 'compact')

nexttile
plot(k_resolutions, Uf_converged, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
xlabel('Number of points for k ($N_k$)', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Critical Flutter Speed $U_f$ (m/s)', 'Interpreter', 'latex', 'FontSize', 12)
title('Convergence vs. Number of Points', 'FontSize', 12)
grid on; grid minor;

nexttile
plot(step_sizes, Uf_converged, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'Color', 'r')
set(gca, 'XDir', 'reverse') % Reverse axis to show increasing resolution left to right
xlabel('Step size $\Delta k$', 'Interpreter', 'latex', 'FontSize', 12)
title('Convergence vs. Step Size', 'FontSize', 12)
grid on; grid minor;

fprintf('\nConvergence Study Complete.\n');
fprintf('The base model (ii=150) computed Uf = %.4f m/s.\n', Uf_converged(k_resolutions == 150));
fprintf('The finest model (ii=1000) computed Uf = %.4f m/s.\n', Uf_converged(end));