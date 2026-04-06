clear all;

% Basic wing planform
cr = 1.5 ;         %chord
lambda = 0.4 ;     % Taper ratio
L = 6 ;            % Semi span
e = 0.25 ;         % Elastic axis location relative to quarter chord
LAM = 25 ;         % sweep in degrees

% Meshing part
NP=100; % Number of panels
DY=2/(NP); % Dimensionless panel span
y_L=0:DY:1; % Panel position along one wing
c_Y=cr*(1-lambda*y_L); % chord along the span due to taper
dy=DY*L; % Actual panel span
S=trapz(c_Y)*dy*2; % Wing area
AR=(L*2)^2/S; % Aspect ratio

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

% incident air conditions
mg=10*1000*9.81; %aircraft weight 10 ton aircraft
Vinf=150; % m/s
rho=0.5238*1.225; % kg/m^3 at 20,000 ft
% Elastic properties
EI=2e6; % Nm^2
GJ=5e5; %Nm^2/rad

a3=2*pi*AR/(2+sqrt(4+AR^2)); % Lift Cyrve slope 3D approximation

% number of bending and torsional model
N=3; %Bending polynomials
M=2; %Torsion polynomials

% FUNCTION FILE CREATED TO CALCULATE STIFFNESS MATRIX
E=stiffness_matrix(N,M,EI,GJ,L,y_L)

% step 1
trimstep=1;
% original alpha to support weight if there is no bend or twist
alpha(trimstep) =mg/(0.5*rho*Vinf^2*S*a3);
% First load estimation
Li=0.5*rho*Vinf^2*(cr*(1-lambda*(y_L)))*a3*alpha(trimstep);
% Rigid wing CL along span
CL_y_R=Li./(0.5*rho*Vinf.^2*(cr*(1-lambda*y_L)));

% step 2 Determine forces on each basis function via virtual work
Mi=e*cr*Li;
for ii=1:N
    psi_i(ii,:) =(y_L).^(ii+1); % ith bending function wing
    psi_id(ii,:)=((ii+1)/L)*((y_L).^ii); % first derivative of the ith bending ufnction of the wiing
    F(ii,:)=trapz(y_L,Li.*(psi_i(ii,:)))*L;
end
for ii=1:M
    phi_i(ii,:)= (y_L).^(ii); % ith torsion function wing
    F(ii+N,:)=trapz(y_L,Mi.*(phi_i(ii,:)))*L;
end

% step 3 determine bend and twist from this load by solving system of equations
eta=E\F;
% spanwise twist and bending for the first trim step
theta=phi_i'*eta(N+1:N+M);
wd=psi_id'*eta(1:N);

% step 4 determine angle of attack including elastic effects
alphae=alpha(trimstep)+theta*cosd(LAM)-wd*sind(LAM);

%step 5 recompure lift along the span
Li=(0.5*rho*Vinf^2*(cr*(1-lambda*(y_L)))*a3).*alphae';
Ltot=trapz(y_L,Li)*L*2; % calculate total wing lift
formatSpec = 'TrimStep %i, Ltot %f, mg %f, Lover %f, alpha %f, alphae(NP) %f\n';
fprintf (formatSpec,trimstep, Ltot, mg, Ltot-mg, alpha, alphae(NP/2))

% step 6 use simple trimming routine
for trimstep=2:200
% elastic angle of attack for one wing
alphae=alpha(trimstep-1)+theta*cosd(LAM)-wd*sind(LAM);
% lift for elastic wing
Li=(0.5*rho*Vinf^2*(cr*(1-lambda*(y_L)))*a3).*alphae';
% total lift
Ltot=trapz(y_L,Li)*L*2;
% compute integral of lift and Lift "error"
Lover=Ltot-mg;
% root angle of attack over/under correction
aover=Lover/(.5*rho*Vinf^2*S*a3);

% step 7 reduce alpha root so that the increase difference Lover is reduced
alpha(trimstep)=alpha(trimstep-1)-aover*0.1;
formatSpec = 'TrimStep %i, Ltot %f, mg %f, Lover %f, aover %f, alpha %f\n';
fprintf (formatSpec,trimstep, Ltot, mg, Lover, aover, alpha(trimstep))

% step 8 iterate until a/c model is trimed
Mi=e*cr*Li;
clear psi_i psi_id phi_i F
for ii=1:N
    psi_i(ii,:) =(y_L).^(ii+1); % ith bending function wing
    psi_id(ii,:)=((ii+1)/L)*((y_L).^ii); % first derivative of the ith bending ufnction of the wiing
    F(ii,:)=trapz(y_L,Li.*(psi_i(ii,:)))*L;
end
for ii=1:M
    phi_i(ii,:)= (y_L).^(ii); % ith torsion function wing
    F(ii+N,:)=trapz(y_L,Mi.*(phi_i(ii,:)))*L;
end
% solve system of equations for deflections "eta"
eta=E\F;
% spanwise bending and twist
theta=phi_i'*eta(N+1:N+M);
wd=psi_id'*eta(1:N);

% if within 1 percent of weight don't trim any more
if abs(Lover)<mg/100, break, end % End loop condition

end

% yotal lift coefficient
CL=Ltot/(0.5*rho*(Vinf^2)*S);
% spanwise lift coefficient
CL_y=Li./(0.5*rho*(Vinf^2)*(cr*(1-lambda*(y_L))));
CL_sweep=(CL_y/CL)';

%% Plotting
close all ;
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
plot(y_L, CL_y, '-o',  'LineWidth', 2, 'DisplayName','CL (elastic VLM)', 'Color', [0 0.427 0.831])
plot(y_L, CL_sweep, '-o',  'LineWidth', 2, 'DisplayName', 'CL/CL_{total} (normalised)', 'Color', [0.286 0.678 0])
plot(y_L, CL_y_R, '--r', 'LineWidth', 1.5, 'DisplayName','CL (rigid strip theory seed)')
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Lift coefficient C_L', 'FontSize', 14)
title(sprintf('Spanwise Lift Coefficient – Elastic Wing (LAM = %.0f)', LAM), 'FontSize', 14)
legend('Location','southwest', 'FontSize', 11)

% Elastic angle of attack
figure('Name','Elastic Angle of Attack')
hold on 
grid on
grid minor
plot(y_L, alphae*180/pi, '-o', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Aeroelastic angle (deg)', 'FontSize', 14)
title('Spanwise Elastic Angle of Attack', 'FontSize', 14)

% Trim convergence 
figure('Name','Trim Convergence')
hold on
grid on
grid minor
plot(alpha*180/pi, '-o', 'LineWidth', 2)
xlabel('Iteration', 'FontSize', 14)
ylabel('Root pitch angle (deg)', 'FontSize', 14)
title('Trim Convergence History', 'FontSize', 14)
xlim([0 max(trimstep)])

% Torsion angle 
figure('Name','Torsion Distribution')
hold on 
grid on
grid minor
plot(y_L, theta*180/pi, '-o', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Torsion angle (deg)', 'FontSize', 14)
title('Spanwise Torsion Distribution', 'FontSize', 14)

% Bending slope 
figure('Name','Bending Slope Distribution')
hold on 
grid on
grid minor
plot(y_L, wd, '-o', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Bending slope dw/dy (m)', 'FontSize', 14)
title('Spanwise Bending Slope Distribution', 'FontSize', 14)

% Spanwise loading Li 
figure('Name','Spanwise Loading (VLM elastic)')
hold on 
grid on
grid minor
plot(y_L, Li, '-o', 'LineWidth', 2)
xlabel('Dimensionless span (y/L)', 'FontSize', 14)
ylabel('Local lift L_i (N)', 'FontSize', 14)
title('Spanwise Load Distribution – Elastic Wing (VLM)', 'FontSize', 14)