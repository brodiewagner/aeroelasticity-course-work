%% clean workspace
clear
clc

%% Initialise Parameters
% Basic wing planform 
c = 1.5 ;               % chord
lambda = 2/3 ;          % Taper ratio
L = sqrt(2) ;           % Semi span
e = 0.25 ;              % Elastic axis location relative to quarter chord
LAM = 45 ;              % sweep in degrees

% Meshing
NP = 50 ;                           % Number of panels
DY = 1/NP ;                         % Dimensionless panel span
y_L = 0:DY:1 ;                      % 1x51 panel EDGES (dimensionless)
c_Y = c*(1 - lambda*y_L) ;          % chord at each edge (1x51)
dy = DY*L ;                         % actual panel span
S = trapz(c_Y)*dy*2 ;               % wing area
AR = (L*2)^2/S ;                    % aspect ratio

% 50-point panel CENTRES (used everywhere inside loops)
y_ctrl = 0.5*(y_L(1:NP) + y_L(2:NP+1)) ;    % 1x50, midpoints
c_ctrl = c*(1 - lambda*y_ctrl) ;            % chord at each centre (1x50)

% Incident air conditions
mg = 10*1000*9.81 ;         % aircraft weight (N)
Vinf = 150 ;                % m/s
rho = 0.5238*1.225 ;        % kg/m^3 at 20,000 ft

% Elastic properties
EI = 2e6 ;      % Nm^2
GJ = 5e5 ;      % Nm^2/rad

% 3D lift curve slope estimate
a3 = 2*pi*AR/(2+sqrt(4+AR^2)) ;

% Number of bending and torsional modes
N = 3 ;     % Bending polynomials
M = 2 ;     % Torsion polynomials

% Stiffness matrix (uses y_L edges internally)
E = stiffness_matrix(N, M, EI, GJ, L, y_L) ;

% Wing geometry for VLM (all 50-panel quantities)
x_hc = y_L*L*tand(LAM) ;                % half chord x at edges (1x51)
x_qc = -0.25*c_Y + x_hc ;               % quarter chord x at edges (1x51)
x_3qc =  0.25*c_Y + x_hc ;              % 3/4 chord x at edges (1x51)
z = zeros(size(x_hc)) ;

% VLM points: A = left corner, B = right corner, C = control point (all 3x50)
A = [x_qc(1:NP)  ; y_L(1:NP)*L  ; z(1:NP)  ] ;
B = [x_qc(2:NP+1); y_L(2:NP+1)*L; z(2:NP+1)] ;
Cx_ctrl = interp1(y_L, x_3qc, y_ctrl) ;    % 1x50
C = [Cx_ctrl ; y_ctrl*L ; zeros(1,NP)] ;   % 3x50

% Unit normals (3x50)
n = zeros(3, NP) ;
for k = 1:NP
    nk = cross(A(:,k)-C(:,k), A(:,k)-B(:,k)) ;
    n(:,k) = nk/norm(nk) ;
end

% Rigid wing CL for reference plot
alpha_rigid = mg/(0.5*rho*Vinf^2*S*a3) ;
L_i_R = 0.5*rho*Vinf^2*c_ctrl*a3*alpha_rigid ;   % 1x50
CL_y_R = L_i_R ./ (0.5*rho*Vinf^2*c_ctrl) ;      % 1x50

%% ---- TRIM LOOP ----

% Step 1: initial alpha guess, initial elastic deflections = 0
alpha(1) = alpha_rigid ;
theta = zeros(NP, 1) ;      % 50x1 twist at control points
wd = zeros(NP, 1) ;         % 50x1 bending slope at control points

for trimstep = 1:200

    % Aeroelastic angle of attack at 50 control points (column vector 50x1)
    alphae = alpha(trimstep) + theta*cosd(LAM) - wd*sind(LAM) ;  % 50x1

    % VLM lift (returns 1x50 row vector)
    [L_i, Cltot] = solve_VLM_2(alphae', Vinf, rho, S, NP, A, B, C, n, DY, L) ;

    % Total lift (both wings)
    L_tot = sum(L_i)*2 ;
    L_over = L_tot - mg ;

    % Trim correction
    aover = L_over/(0.5*rho*Vinf^2*S*a3) ;
    
    fprintf('TrimStep %i, Ltot %.2f, mg %.2f, Lover %.2f, alpha %.6f\n', ...
            trimstep, L_tot, mg, L_over, alpha(trimstep)) ;

    if trimstep < 200
        alpha(trimstep+1) = alpha(trimstep) - aover*0.05 ;
    end

    % Virtual work: forces on each basis function (using y_ctrl, 50 points)
    Mi = e*c_ctrl.*L_i ;    % moment (1x50)

    % Build basis function matrices (N+M) x 50
    psi_i = zeros(N, NP) ;
    psi_id = zeros(N, NP) ;
    phi_i = zeros(M, NP) ;
    F = zeros(N+M, 1) ;

    for ii = 1:N
        psi_i(ii,:) = y_ctrl.^(ii+1) ;
        psi_id(ii,:) = ((ii+1)/L)*y_ctrl.^ii ;
        F(ii) = trapz(y_ctrl, L_i.*psi_i(ii,:))*L ;
    end
    for ii = 1:M
        phi_i(ii,:) = y_ctrl.^ii ;
        F(N+ii) = trapz(y_ctrl, Mi.*phi_i(ii,:))*L ;
    end

    % Solve for deflection coefficients
    eta = E\F ;

    % Spanwise bending slope and twist at 50 control points (50x1 columns)
    theta = (phi_i')*eta(N+1:N+M) ;             % 50x1
    wd = (psi_id')*eta(1:N) ;                   % 50x1

    % Convergence check
    if abs(L_over) < mg/100
        fprintf('Converged at step %i\n', trimstep) ;
        break
    end

end

%% ---- POST PROCESSING ----

CL = L_tot/(0.5*rho*Vinf^2*S) ;

% Spanwise CL (1x50)
CL_y = L_i./(0.5*rho*Vinf^2*c_ctrl) ;
CL_sweep = CL_y/CL ;

figure('name', 'Lift Coefficient')
plot(y_ctrl, CL_y, '-b', 'LineWidth', 2) ; 
hold on
plot(y_ctrl, CL_sweep, '-g', 'LineWidth', 2) ;
plot(y_ctrl, CL_y_R, '-r', 'LineWidth', 2) ;
xlabel('Dimensionless span (y/L)') ;
ylabel('Lift coefficient (CL)') ;
legend('CL elastic','CL/CL_{avg}','CL rigid') ;
grid on

figure('Name', 'Aeroelastic Angle')
plot(y_ctrl, alphae*180/pi, '-b', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)')
ylabel('Aeroelastic angle (deg)')
grid on

figure('Name', 'Pitch Angle')
plot(alpha*180/pi, '-b', 'LineWidth', 2)
xlabel('Iteration')
ylabel('Pitch angle (deg)')
grid on

figure('name', 'Torsion Angle')
plot(y_ctrl, theta*180/pi, '-b', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)')
ylabel('Torsion angle (deg)')
grid on

figure('Name', 'Bending Slope')
plot(y_ctrl, wd, '-b', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)')
ylabel('Bending slope (m/m)')
grid on