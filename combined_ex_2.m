% ==========================================================================================================================
% Aeroelastic trim of elastic wing using Vortex Lattice Method (VLM) 
% for aerodynamic loads instead of strip theory.
%
% The VLM solves the full 3-D induced downwash problem across all panels 
% simultaneously, so spanwise load redistribution due to sweep, taper and 
% elastic twist/bending is captured correctly.
%
% MESH CONVENTION (important):
%   solve_VLM / vortex_lattice_method_rigid use a FULL-SPAN mesh with NP
%   panels running tip-to-tip (y/L in [-1, +1]).  The structural solver
%   works on ONE wing (y/L in [0, 1]) with NP/2 + 1 stations.
%   We therefore:
%     - Build a full-span VLM mesh with NP_vlm = NP_struct * 2 panels.
%     - Mirror the elastic alpha distribution symmetrically onto both 
%       wings before each VLM call.
%     - Extract only the starboard half of the VLM lift result for the
%       structural integrals (panels NP_vlm/2+1 : NP_vlm).
%
% REQUIREMENTS:
%   V_AB.m, VA_INF.m, VB_INF.m, solve_VLM.m, stiffness_matrix.m must 
%   all be on the MATLAB path.
%
% PROBLEMS:
%   Values for Lift Coefficient are a little low at the root. dk why
% ==========================================================================================================================

clear
clc

%% WING PLANFORM PARAMETERS
cr = 1.5 ;          % root chord (m)
lambda = 0.4 ;      % taper ratio, ct/cr
L = 6 ;             % semi-span (m)
e = 0.25 ;          % elastic axis offset aft of quarter-chord (x chords)
LAM = 0 ;           % half-chord sweep angle (degrees)

%% STRUCTURAL MESH (half-span, NP_s stations)
NP_s = 50 ;                         % half-span structural panels (keep even)
DY_s = 1/NP_s ;                     % dimensionless panel span (one wing)
y_L = 0:DY_s:1 ;                    % NP_s+1 stations, y/L in [0,1]
c_Y = cr*(1 - lambda*y_L) ;         % chord distribution along half-span
dy = DY_s*L ;                       % physical panel span (m)
S = trapz(c_Y)*dy*2 ;               % wing reference area (trapz, both wings)
AR = (2*L)^2/S ;                    % aspect ratio

%%  VLM MESH (full-span, NP_v panels)
% The VLM functions expect NP panels spanning tip-to-tip with y/L running from -1 to +1. NP_v MUST be even.
NP_v = 2*NP_s ;             % total full-span VLM panels
DY_v = 1/NP_s ;             % dimensionless VLM panel span
Y_L_v = 0:DY_v:1 ;          % half-span ordinate for VLM chord calc
y_L_v = -1:DY_v:1 ;         % full-span ordinate (NP_v+1 points)

% Chord length across span, tip-to-tip (needed for x-coordinates)
C_Y_v = cr*(1 - lambda*Y_L_v) ;         % for half-span
c_Y_v = [C_Y_v(end:-1:2), C_Y_v] ;      % mirror: left tip -> right tip

% Geometry lines (half-chord sweep reference)
x_hc = abs(y_L_v) * L * tand(LAM);                  % half-chord x
x_qc = -0.25*c_Y_v + x_hc ;                         % quarter-chord x
x_le = -0.5*c_Y_v+x_hc ;                            % leading edge
x_te = 0.5*c_Y_v+x_hc ;                             % trailing edge
x_3qc =  0.25*c_Y_v + x_hc ;                        % 3/4-chord x
z = zeros(size(y_L_v)) ;

% Vortex corners (A = left end, B = right end of bound vortex)
A = [x_qc(1:NP_v); y_L_v(1:NP_v)*L; z(1:NP_v)] ;
B = [x_qc(2:NP_v + 1); y_L_v(2:NP_v + 1)*L; z(2:NP_v + 1)] ;

% Collocation points at 3/4-chord, mid-panel spanwise
Cy = 0.5*(A(2,:) + B(2,:)) ;
Cx = interp1(y_L_v, x_3qc, Cy/L) ;
C = [Cx; Cy; 0*Cy] ;

% Panel unit normals (flat wing -> purely z-direction, but computed properly so code generalises to cambered/dihedral cases)
n_v = zeros(3, NP_v) ;
for k = 1:NP_v
    nk = cross(A(:,k) - C(:,k), A(:,k) - B(:,k)) ;
    n_v(:,k) = nk / norm(nk) ;
end

fprintf('VLM mesh built: %d full-span panels  (NP_s=%d half-span structural panels)\n', NP_v, NP_s);

%%  FLIGHT CONDITIONS
mg = 10e3*9.81 ;            % aircraft weight (N)  
Vinf = 150 ;                % airspeed (m/s)
rho = 0.5238*1.225 ;        % air density at ~20 000 ft (kg/m^3)
qinf = 0.5*rho*Vinf^2 ;     % dynamic pressure (Pa)

%%  STRUCTURAL MODEL
EI = 2e6 ;          % bending stiffness (N·m^2)
GJ = 5e5 ;          % torsional stiffness (N·m^2/rad)
N = 3 ;             % number of bending polynomial modes
M = 2 ;             % number of torsion polynomial modes

% Structural stiffness matrix 
E_stiff = stiffness_matrix(N, M, EI, GJ, L, y_L) ;

tic
%%  INITIAL TRIM ESTIMATE (use strip theory once for step 1 seed)
% 3-D lift curve slope (Helmbold formula) used only for first guess of root alpha so the trim loop starts near the solution.
a3 = 2*pi*AR / (2 + sqrt(4 + AR^2)) ;

trimstep = 1 ;
alpha(trimstep) = mg/(qinf*S*a3) ;                  % rigid-wing root alpha estimate

% Rigid-wing lift (strip theory seed, half-span) 
Li_s = qinf*c_Y*a3*alpha(trimstep) ;                % 1 x (NP_s+1) half-span lift/m * chord

% Store rigid CL distribution for final comparison plot
CL_y_R = Li_s./(qinf*c_Y) ;                         % local CL, rigid strip theory

%% build basis functions on the structural half-span grid
% declare them outside the loop 
psi_i = zeros(N, NP_s+1) ;
psi_id = zeros(N, NP_s+1) ;
phi_i = zeros(M, NP_s+1) ;
Mi_s = e*cr*Li_s ;            % pitching moment per unit span (strip seed)
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

%% STEP 3 – Structural solve
eta = E_stiff\F ;
theta = phi_i'*eta(N+1:N+M) ;           % (NP_s+1) x 1 torsion
wd = psi_id'*eta(1:N) ;                 % (NP_s+1) x 1 bending slope

%% STEP 4 – Elastic angle of attack (half-span)
alphae_half = alpha(trimstep) + theta*cosd(LAM) - wd*sind(LAM) ;        % (NP_s+1) x 1

%%  STEP 5 – First VLM call with elastic alpha
% Mirror alphae onto full span: port wing is symmetric
%   Full-span panels 1..NP_v/2 = port(y < 0), mirrored from starboard
%   Full-span panels NP_v/2+1..NP_v = starboard (y > 0)

% Panel centres for the structural grid (NP_s values, one per panel)
alphae_panels = 0.5*(alphae_half(1:end-1) + alphae_half(2:end)) ;       % NP_s x 1

% Mirror to full span (port = flip of starboard)
alphae_full = [flip(alphae_panels); alphae_panels]' ;                   % 1 x NP_v  row vector
 
% VLM call
[Li_vlm, Cltot] = solve_VLM(alphae_full, Vinf, rho, S, NP_v, A, B, C, n_v, DY_v, L) ;
Li_vlm = Li_vlm(:)' ;                                                   % force 1 x NP row vector

% Extract starboard half-span lift (panels NP_v/2+1 : NP_v)
Li_star = Li_vlm(NP_v/2+1 : NP_v)/dy ;                                     % 1 x NP_s  (one value per structural panel)

% Map panel lift onto structural stations by linear interpolation (structural stations are at panel edges; panels are centred between them)
y_panels = 0.5*(y_L(1:end-1) + y_L(2:end)) ;                            % NP_s panel-centre y/L
Li_s = interp1(y_panels, Li_star, y_L, 'linear') ;                      % 1 x (NP_s+1)

L_tot = sum(Li_vlm) ;
fprintf('Step %3i | L_tot = %8.1f N | mg = %12.1f N | error = %.1f N | alpha = %.4f degrees\n', ...
        trimstep, L_tot, mg, L_tot-mg, alpha(trimstep)*180/pi) ;

%%  TRIM LOOP (Steps 6–8)
for trimstep = 2:300

    % STEP 6: update elastic AoA from previous structural state 
    alphae_half = alpha(trimstep-1) + theta*cosd(LAM) - wd*sind(LAM) ;                  % (NP_s+1) x 1

    % Mirror to full span (row vector, NP_v values)
    alphae_panels = 0.5*(alphae_half(1:end-1) + alphae_half(2:end)) ;
    alphae_full   = [flip(alphae_panels) ; alphae_panels]' ;

    % VLM solve (3-D, accounts for sweep/taper/twist effects) 
    [Li_vlm, Cltot] = solve_VLM(alphae_full, Vinf, rho, S, NP_v, A, B, C, n_v, DY_v, L) ;

    % Starboard half only
    Li_star = Li_vlm(NP_v/2+1 : NP_v)/dy ;

    % Interpolate panel lift to structural stations
    Li_s = interp1(y_panels, Li_star, y_L, 'linear', 'extrap') ;

    % Total lift
    L_tot  = sum(Li_vlm) ;
    L_over = L_tot - mg ;

    % STEP 7: trim correction to root alpha 
    % Correction proportional to lift error; scale by full-span VLM CL_alpha
    % (use qinf*S*Cltot/alpha as an approximate local sensitivity)
    a_over = L_over / (qinf * S * a3) ;
    alpha(trimstep) = alpha(trimstep-1) - a_over * 0.1 ;

    fprintf('Step %3i | L_tot = %8.1f N | L_over = %+8.1f N | a_over = %+.5f | alpha = %.5f degrees\n', ...
            trimstep, L_tot, L_over, a_over, alpha(trimstep)*180/pi) ;

    % STEP 8: structural solve with updated aerodynamic loads 
    Mi_s = e * cr * Li_s ;

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

    % Convergence check: within 1% of weight
    if abs(L_over) < mg*0.01
        fprintf('\nConverged at trim step %i\n', trimstep) ;
        break
    end
end

%%  POST-PROCESSING
% Final elastic AoA
alphae_half = alpha(end) + theta*cosd(LAM) - wd*sind(LAM) ;

% Final spanwise lift coefficient
CL_y = Li_s./(qinf*c_Y) ;               % local CL
CL_total = L_tot/(qinf*S) ;              % wing CL
CL_norm = (CL_y/CL_total)' ;            % normalised

toc
%%  PLOTS

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

% Spanwise lift coefficient 
figure('Name','Spanwise Lift Coefficient (VLM elastic)')
hold on 
grid on
grid minor
plot(y_L, CL_y, '-b',  'LineWidth', 2, 'DisplayName','CL (elastic VLM)')
plot(y_L, CL_norm, '-g',  'LineWidth', 2, 'DisplayName','CL/CL_{total} (normalised)')
plot(y_L, CL_y_R, '-r', 'LineWidth', 1.5, 'DisplayName','CL (rigid strip theory seed)')
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Lift coefficient C_L', 'FontSize', 14)
title(sprintf('Spanwise Lift Coefficient – Elastic Wing (LAM = %.0f)', LAM), 'FontSize', 14)
legend('Location','southwest', 'FontSize', 11)

% % Elastic angle of attack
% figure('Name','Elastic Angle of Attack')
% hold on 
% grid on
% grid minor
% plot(y_L, alphae_half*180/pi, '-b', 'LineWidth', 2)
% xlabel('Dimensionless span (y/L)', 'FontSize', 14)
% ylabel('Aeroelastic angle (deg)', 'FontSize', 14)
% title('Spanwise Elastic Angle of Attack', 'FontSize', 14)
% 
% % Trim convergence 
% figure('Name','Trim Convergence')
% hold on
% grid on
% grid minor
% plot(alpha*180/pi, '-b', 'LineWidth', 2)
% xlabel('Iteration', 'FontSize', 14)
% ylabel('Root pitch angle (deg)', 'FontSize', 14)
% title('Trim Convergence History', 'FontSize', 14)
% xlim([0 max(trimstep)])
% 
% % Torsion angle 
% figure('Name','Torsion Distribution')
% hold on 
% grid on
% grid minor
% plot(y_L, theta*180/pi, '-b', 'LineWidth', 2)
% xlabel('Dimensionless span (y/L)', 'FontSize', 14)
% ylabel('Torsion angle (deg)', 'FontSize', 14)
% title('Spanwise Torsion Distribution', 'FontSize', 14)
% 
% % Bending slope 
% figure('Name','Bending Slope Distribution')
% hold on 
% grid on
% grid minor
% plot(y_L, wd, '-b', 'LineWidth', 2)
% xlabel('Dimensionless span (y/L)', 'FontSize', 14)
% ylabel('Bending slope dw/dy (m)', 'FontSize', 14)
% title('Spanwise Bending Slope Distribution', 'FontSize', 14)
% 
% % Spanwise loading Li 
% figure('Name','Spanwise Loading (VLM elastic)')
% hold on 
% grid on
% grid minor
% plot(y_L, Li_s, '-b', 'LineWidth', 2)
% xlabel('Dimensionless span (y/L)', 'FontSize', 14)
% ylabel('Local lift L_i (N)', 'FontSize', 14)
% title('Spanwise Load Distribution – Elastic Wing (VLM)', 'FontSize', 14)