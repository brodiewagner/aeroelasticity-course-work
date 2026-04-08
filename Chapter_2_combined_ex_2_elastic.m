clear
clc
% close all

%% initialise variables
% wing parameters
cr = 1.5 ;          % root chord (m)
lambda = 2/3 ;      % taper ratio  ct/cr
L = 6 ;             % semi-span (m)
e = 0.25 ;          % elastic axis offset aft of quarter-chord (x chords)
LAM = -45 ;           % half-chord sweep angle (degrees)

% flight conditions
mg = 10e3*9.81 ;            % aircraft weight (N)  
Vinf = 150 ;                % airspeed (m/s)
rho = 0.5238*1.225 ;        % air density at ~20 000 ft (kg/m^3)
qinf = 0.5*rho*Vinf^2 ;     % dynamic pressure (Pa)

% structural parameters
EI = 2e6 ;          % bending stiffness (N·m^2)
GJ = 5e5 ;          % torsional stiffness (N·m^2/rad)

% polynomials
N = 3 ;             % number of bending polynomial modes
M = 2 ;             % number of torsion polynomial modes

%% meshing 
% structural mesh (half-span, NP_s stations)
NP_s = 50 ;                         % half-span structural panels (keep even)
DY_s = 1/NP_s ;                     % dimensionless panel span (one wing)
y_L = 0:DY_s:1 ;                    % NP_s+1 stations, y/L in [0,1]
c_Y = cr*(1 - lambda*y_L) ;         % chord distribution along half-span
dy = DY_s*L ;                       % physical panel span (m)
S = trapz(c_Y)*dy*2 ;               % wing reference area (trapz, both wings)
AR = (2*L)^2/S ;                    % aspect ratio

% structural stiffness matrix 
E_stiff = stiffness_matrix(N, M, EI, GJ, L, y_L) ;

% VLM mesh 
NP_v = 2*NP_s ;             % total full-span VLM panels
DY_v = 1/NP_s ;             % dimensionless VLM panel span
Y_L_v = 0:DY_v:1 ;          % half-span ordinate for VLM chord calc
y_L_v = -1:DY_v:1 ;         % full-span ordinate (NP_v+1 points)

% chord length across span, tip-to-tip
C_Y_v = cr*(1 - lambda*Y_L_v) ;         % for half-span
c_Y_v = [C_Y_v(end:-1:2), C_Y_v] ;      % mirror: left tip -> right tip

% geometry lines 
x_hc = abs(y_L_v) * L * tand(LAM);                  % half-chord x
x_qc = -0.25*c_Y_v + x_hc ;                         % quarter-chord x
x_le = -0.5*c_Y_v+x_hc ;                            % leading edge
x_te = 0.5*c_Y_v+x_hc ;                             % trailing edge
x_3qc =  0.25*c_Y_v + x_hc ;                        % 3/4-chord x
z = zeros(size(y_L_v)) ;

% vortex corners (A = left end, B = right end of bound vortex)
A = [x_qc(1:NP_v); y_L_v(1:NP_v)*L; z(1:NP_v)] ;
B = [x_qc(2:NP_v + 1); y_L_v(2:NP_v + 1)*L; z(2:NP_v + 1)] ;

% collocation points at 3/4-chord, mid-panel spanwise
Cy = 0.5*(A(2,:) + B(2,:)) ;
Cx = interp1(y_L_v, x_3qc, Cy/L) ;
C = [Cx; Cy; 0*Cy] ;

% panel unit normals 
n_v = zeros(3, NP_v) ;
for k = 1:NP_v
    nk = cross(A(:,k) - C(:,k), A(:,k) - B(:,k)) ;
    n_v(:,k) = nk / norm(nk) ;
end

tic

%% intial rigid trim 
% intial trim estimate 
a3 = 2*pi*AR / (2 + sqrt(4 + AR^2)) ;

trimstep = 1 ;
alpha(trimstep) = mg/(qinf*S*a3) ;                  % rigid-wing root alpha estimate

% rigid-wing lift (strip theory seed, half-span) 
Li_s = qinf*c_Y*a3*alpha(trimstep) ;                % 1 x (NP_s+1) half-span lift/m * chord

% store rigid CL distribution for final comparison plot
CL_y_R = Li_s./(qinf*c_Y) ;                         % local CL, rigid strip theory

% basis functions on the structural half-span grid
psi_i = zeros(N, NP_s+1) ;
psi_id = zeros(N, NP_s+1) ;
phi_i = zeros(M, NP_s+1) ;
Mi_s = e.*c_Y.*Li_s ;            % pitching moment per unit span (strip seed)
F = zeros(N+M, 1) ;

for ii = 1:N
    psi_i(ii,:) = y_L.^(ii+1) ;
    psi_id(ii,:) = ((ii+1)/L)*y_L.^ii ;
    F(ii) = trapz(y_L, Li_s.*psi_i(ii,:))*L ;
end

for ii = 1:M
    phi_i(ii,:) = y_L.^ii ;
    F(ii+N) = trapz(y_L, Mi_s.*phi_i(ii,:))*L ;
end

% structural solver
eta = E_stiff\F ;
theta = phi_i'*eta(N+1:N+M) ;           % (NP_s+1) x 1 torsion
wd = psi_id'*eta(1:N) ;                 % (NP_s+1) x 1 bending slope

% elastic angle of attack (half-span)
alphae_half = alpha(trimstep) + theta*cosd(LAM) - wd*sind(LAM) ;        % (NP_s+1) x 1

%% first VLM call with elastic alpha
% mirror alphae onto full span: port wing is symmetric

% panel centres for the structural grid (NP_s values, one per panel)
alphae_panels = 0.5*(alphae_half(1:end-1) + alphae_half(2:end)) ;       % NP_s x 1

% mirror to full span (port = flip of starboard)
alphae_full = [flip(alphae_panels); alphae_panels]' ;                   % 1 x NP_v  row vector
 
% VLM call
[Li_vlm] = solve_VLM(alphae_full, Vinf, rho, S, NP_v, A, B, C, n_v, DY_v, L) ;
Li_vlm = Li_vlm(:)' ;                                                   % force 1 x NP row vector

% extract starboard half-span lift (panels NP_v/2+1 : NP_v)
Li_star = Li_vlm(NP_v/2+1 : NP_v)/dy ;                                     % 1 x NP_s  (one value per structural panel)

% map panel lift onto structural stations by linear interpolation (structural stations are at panel edges; panels are centred between them)
y_panels = 0.5*(y_L(1:end-1) + y_L(2:end)) ;                            % NP_s panel-centre y/L
Li_s = interp1(y_panels, Li_star, y_L, 'linear') ;                      % 1 x (NP_s+1)

L_tot = sum(Li_vlm) ;
fprintf('Step %3i | L_tot = %8.1f N | mg = %12.1f N | error = %.1f N | alpha = %.4f degrees\n', ...
        trimstep, L_tot, mg, L_tot-mg, alpha(trimstep)*180/pi) ;

%%  trim loop
for trimstep = 2:300

    % update elastic AoA from previous structural state 
    alphae_half = alpha(trimstep-1) + theta*cosd(LAM) - wd*sind(LAM) ;                  % (NP_s+1) x 1

    % mirror to full span (row vector, NP_v values)
    alphae_panels = 0.5*(alphae_half(1:end-1) + alphae_half(2:end)) ;
    alphae_full   = [flip(alphae_panels) ; alphae_panels]' ;

    % VLM solve 
    [Li_vlm] = solve_VLM(alphae_full, Vinf, rho, S, NP_v, A, B, C, n_v, DY_v, L) ;

    % starboard half only
    Li_star = Li_vlm(NP_v/2+1 : NP_v)/dy ;

    % interpolate panel lift to structural stations
    Li_s = interp1(y_panels, Li_star, y_L, 'linear', 'extrap') ;

    % total lift
    L_tot  = sum(Li_vlm) ;
    L_over = L_tot - mg ;

    % trim correction to root alpha 
    a_over = L_over / (qinf * S * a3) ;
    alpha(trimstep) = alpha(trimstep-1) - a_over * 0.1 ;

    fprintf('Step %3i | L_tot = %8.1f N | L_over = %+8.1f N | a_over = %+.5f | alpha = %.5f degrees\n', ...
            trimstep, L_tot, L_over, a_over, alpha(trimstep)*180/pi) ;

    % structural solve with updated aerodynamic loads 
    Mi_s = e .* c_Y .* Li_s ;

    F = zeros(N+M, 1) ;
    for ii = 1:N
        F(ii)   = trapz(y_L, Li_s .* psi_i(ii,:))  * L ;
    end
    for ii = 1:M
        F(ii+N) = trapz(y_L, Mi_s .* phi_i(ii,:)) * L ;
    end

    eta   = E_stiff \ F ;
    theta = phi_i' * eta(N+1:N+M) ;
    wd    = psi_id' * eta(1:N)     ;

    % convergence check: within 1% of weight
    if abs(L_over) < mg*0.01
        fprintf('\nConverged at trim step %i\n', trimstep) ;
        break
    end
end

%% post process
% final elastic AoA
alphae_half = alpha(end) + theta*cosd(LAM) - wd*sind(LAM) ;

% final spanwise lift coefficient
CL_y = Li_s./(qinf*c_Y) ;               % local CL
CL_total = L_tot/(qinf*S) ;              % wing CL
CL_norm = (CL_y/CL_total)' ;            % normalised

toc

%% plots
% wing geometry plotting
figure('Name', 'Wing')
title(sprintf('Wing Profile (LAM = %.0f)', LAM), 'FontSize', 14)
plot3(x_le,y_L_v*L,z) ;                                 % plot leading edge
hold on
plot3(x_hc,y_L_v*L,z,'--') ;                            % plot half chord
plot3(x_te,y_L_v*L,z) ;                                 % plot trailing edge
plot3([x_le(NP_v+1) x_te(NP_v+1)],[L L],[0 0]) ;        % plot chord at tips
plot3([x_le(1) x_te(1)],[-L -L],[0 0]) ;                % plot chord at tips
plot3(x_qc,y_L_v*L,z,'r--')                             % plot out quarter chord, along which vortex positions lie
plot3(x_3qc,y_L_v*L,z,'g--')                            % plot out 3/4 chord, along which collocation points lie
grid on 
axis equal
% plot points A,B,C as a check
plot3(A(1,:), A(2,:), A(3,:),'ro')
plot3(B(1,:), B(2,:), B(3,:),'r+')
plot3(C(1,:), C(2,:), C(3,:),'go')
hold off

% spanwise lift coefficient 
figure('Name','Spanwise Lift Coefficient (VLM elastic)')
hold on 
grid on
grid minor
plot(y_L, CL_y, '-o',  'LineWidth', 2, 'DisplayName','CL (elastic VLM)', 'Color', [0 0.427 0.831])
plot(y_L, CL_norm, '-o',  'LineWidth', 2, 'DisplayName', 'CL/CL_{total} (normalised)', 'Color', [0.286 0.678 0])
plot(y_L, CL_y_R, '--r', 'LineWidth', 1.5, 'DisplayName','CL (rigid strip theory seed)')
set(gca, 'fontsize', 20)
xlabel('Dimensionless span (y/L)', 'FontSize', 25)
ylabel('Lift coefficient C_L', 'FontSize', 25)
title(sprintf('Spanwise Lift Coefficient – Elastic Wing (LAM = %.0f)', LAM), 'FontSize', 20)
legend('Location','southwest', 'FontSize', 25)


% elastic angle of attack
figure('Name','Elastic Angle of Attack')
hold on 
grid on
grid minor
plot(y_L, alphae_half*180/pi, '-o', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Aeroelastic angle (deg)', 'FontSize', 14)
title('Spanwise Elastic Angle of Attack', 'FontSize', 14)

% trim convergence 
figure('Name','Trim Convergence')
hold on
grid on
grid minor
plot(alpha*180/pi, '-o', 'LineWidth', 2)
xlabel('Iteration', 'FontSize', 14)
ylabel('Root pitch angle (deg)', 'FontSize', 14)
title('Trim Convergence History', 'FontSize', 14)
xlim([0 max(trimstep)])

% torsion angle 
figure('Name','Torsion Distribution')
hold on 
grid on
grid minor
plot(y_L, theta*180/pi, '-o', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Torsion angle (deg)', 'FontSize', 14)
title('Spanwise Torsion Distribution', 'FontSize', 14)

% bending slope 
figure('Name','Bending Slope Distribution')
hold on 
grid on
grid minor
plot(y_L, wd, '-o', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Bending slope dw/dy (m)', 'FontSize', 14)
title('Spanwise Bending Slope Distribution', 'FontSize', 14)

% spanwise loading Li 
figure('Name','Spanwise Loading (VLM elastic)')
hold on 
grid on
grid minor
plot(y_L, (Li_s*(10^-3)), '-o', 'LineWidth', 2)
set(gca, 'fontsize', 20)
xlabel('Dimensionless span (y/L)', 'FontSize', 25)
ylabel('Local lift L_i (kN)', 'FontSize', 25)
title('Spanwise Load Distribution – Elastic Wing (VLM)', 'FontSize', 20)
legend('Unswept', 'Aft Sweep', 'Fore Sweep')