clear
clc

% Basic wing planform
c = 1.5 ;           %chord
lambda = 2/3 ;      % Taper ratio
L = sqrt(2) ;       % Semi span
e = 0.25 ;          % Elastic axis location relative to quarter chord
LAM = 45 ;          % sweep in degrees

% Meshing part
NP = 50 ;                       % Number of panels
DY = 1/(NP) ;                   % Dimensionless panel span
y_L = -1:DY:1 ;                 % Panel position along both wings
c_Y = c*(1-lambda*y_L) ;        % chord along the span due to taper
dy = DY*L ;                     % Actual panel span
S = trapz(c_Y)*dy*2 ;           % Wing area
AR = (L*2)^2/S ;                % Aspect ratio

% Incident air conditions
mg = 10*1000*9.81 ;         % aircraft weight 10 ton aircraft
Vinf = 150 ;                % m/s
rho = 0.5238*1.225 ;        % kg/m^3 at 20,000 ft

% Elastic properties
EI = 2e6 ;          % Nm^2
GJ = 5e5 ;          % Nm^2/rad

%m0 = 2*pi ;
a3 = 2*pi*AR/(2+sqrt(4+AR^2)) ;         % Estimate Lift Curve slope 3D approximation (Replace with VLM)
%a3 = m0*cosd(LAM)/(sqrt(1 + (m0*cosd(LAM)/(pi*AR)^2) + (m0*cosd(LAM)/(pi*AR)))) ;         % corrected a3 to consider sweep 

% number of bending and torsional model
N = 3 ;         % Bending polynomials
M = 2 ;         % Torsion polynomials

% FUNCTION FILE CREATED TO CALCULATE STIFFNESS MATRIX
E = stiffness_matrix(N,M,EI,GJ,L,y_L) ;

% plot leading edge, trailing edge, 1/4 chord, 1/2 chord and 3/4 chord axes first half chord axis line from tip to tip:
x_hc = sign(y_L).*y_L*L*tand(LAM) ;             % half chord
x_le = -0.5*c_Y+x_hc ;                          % leading edge
x_te = 0.5*c_Y+x_hc ;                           % trailing edge
x_qc = -0.25*c_Y+x_hc ;                         % quarter chord
x_3qc = 0.25*c_Y+x_hc ;                         % 3/4 chord
z = zeros(size(x_hc)) ;                         % all points in a plane

% compute points "A,B,C" for wing
A = [x_qc(1:NP); y_L(1:NP)*L; z(1:NP)] ;        % left vortex corner
B = [x_qc(2:NP+1); y_L(2:NP+1)*L; z(2:NP+1)] ;  % right vortex corner
Cy = 0.5*(A(2,:)+B(2,:)) ;                      % control points lie inbetween A and B for y coord
Cx = interp1(y_L,x_3qc,Cy/L) ;                  % interpolate x coordinate at 3/4 chord at these y points
C = [Cx; Cy; 0*Cy] ;

% get unit normals for each panel
for k = 1:NP
    n(:,k) = cross(A(:,k)-C(:,k),A(:,k)-B(:,k)) ;   % unit normal from plane formed by A,B,C
    n(:,k) = n(:,k)/sqrt(dot(n(:,k),n(:,k))) ;
end
 
%% step 1
trimstep = 1 ;

    % original alpha to support weight if there is no bend or twist
    alpha(trimstep) = mg/(0.5*rho*Vinf^2*S*a3) ;
    
    % First load estimation (needs to be corrected, shape is wrong)
    L_i = 0.5*rho*Vinf^2*(c*(1-lambda*(y_L)))*a3*alpha(trimstep) ;
    
    % Rigid wing CL along span
    CL_y_R = L_i./(0.5*rho*Vinf.^2*(c*(1-lambda*y_L))) ;

% step 2 Determine forces on each basis function via virtual work
Mi = e*c*L_i ;
for ii = 1:N
    psi_i(ii,:) = (y_L).^(ii+1);                     % ith bending function wing
    psi_id(ii,:) = ((ii+1)/L)*((y_L).^ii);           % first derivative of the ith bending ufnction of the wiing
    F(ii,:) = trapz(y_L,L_i.*(psi_i(ii,:)))*L;
end
for ii = 1:M
    phi_i(ii,:) = (y_L).^(ii) ;                      % ith torsion function wing
    F(ii+N,:) = trapz(y_L,Mi.*(phi_i(ii,:)))*L;
end

% step 3 determine bend and twist from this load by solving system of equations
eta = E\F ;

% spanwise twist and bending for the first trim step
theta = phi_i'*eta(N+1:N+M) ;
wd = psi_id'*eta(1:N) ;

% step 4 determine angle of attack including elastic effects
alphae = alpha(trimstep)+theta*cosd(LAM)-wd*sind(LAM) ;             % e stands for aeroelasticity

%step 5 recompute lift along the span
L_i = (0.5*rho*Vinf^2*(c*(1-lambda*(y_L)))*a3).*alphae' ;
L_tot = trapz(y_L,L_i)*L*2 ;                                          % calculate total wing lift
formatSpec = 'TrimStep %i, Ltot %f, mg %f, L_over %f, alpha %f, alphae(NP) %f\n';
fprintf(formatSpec,trimstep, L_tot, mg, L_tot-mg, alpha, alphae(NP/2))

% step 6 use simple trimming routine
for trimstep = 2:200

    % elastic angle of attack for one wing
    alphae = alpha(trimstep - 1) + theta*cosd(LAM) - wd*sind(LAM) ;
    
    % lift for elastic wing
    alphae_ctrl = alphae(1:NP) ; 
    [Li, Cltot] = solve_VLM(alphae', Vinf, rho, S, NP, A, B, C, n, DY, L) ;      % can use different function (strip theory not fully accurate)
    L_tot = sum(L_i)*2 ;                                                        % total lift
    L_over = L_tot - mg ;                                                       % compute integral of lift and Lift "error"
    
    % root angle of attack over/under correction
    aover = L_over/(0.5*rho*Vinf^2*S*a3) ;
    
    % step 7 reduce alpha root so that the increase difference Lover is reduced
    alpha(trimstep) = alpha(trimstep-1)-aover(1)*0.1 ;
    1 ; %-aover*0.9 ;
    formatSpec = 'TrimStep %i, Ltot %f, mg %f, Lover %f, aover %f, alpha %f\n';
    fprintf (formatSpec,trimstep, L_tot, mg, L_over, aover, alpha(trimstep))
    
    % step 8 iterate until a/c model is trimed
    Mi = e*c*L_i ;

    clear psi_i psi_id phi_i F
    for ii = 1:N
        psi_i(ii,:) = (y_L).^(ii+1) ;                   % ith bending function wing
        psi_id(ii,:) = ((ii+1)/L)*((y_L).^ii) ;         % first derivative of the ith bending ufnction of the wiing
        F(ii,:) = trapz(y_L, L_i.*(psi_i(ii,:)))*L ;
    end
    for ii = 1:M
        phi_i(ii,:) = (y_L).^(ii) ;                     % ith torsion function wing
        F(ii+N,:) = trapz(y_L, Mi.*(phi_i(ii,:)))*L ;
    end
    
    % solve system of equations for deflections "eta"
    eta = E\F ;
    
    % spanwise bending and twist
    theta = phi_i'.*eta(N+1:N+M) ;
    wd = psi_id'.*eta(1:N) ;
    
    % if within 1 percent of weight don't trim any more
    if abs(L_over) < mg/100 
        break 
    end              % End loop condition

end

% yotal lift coefficient
CL = L_tot/(0.5*rho*(Vinf^2)*S) ;

% spanwise lift coefficient
CL_y = L_i./(0.5*rho*(Vinf^2)*(c*(1-lambda*(y_L)))) ;
CL_sweep = (CL_y/CL)' ;

figure
hold on
plot (y_L,CL_y, '-b', 'LineWidth', 2)
hold on
plot(y_L,CL_sweep, '-g', 'LineWidth', 2)
hold on
plot (y_L,CL_y_R(1:51), '-r', 'LineWidth', 2)
xlabel ('Dimensionless span (y/L)')
ylabel ('Lift coefficeint (CL)')
grid on
hold off

figure
hold on
plot (y_L,alphae*180/pi, '-b', 'LineWidth', 2)
xlabel ('Dimensionless span (y/L)')
ylabel ('Aeroelastic angle (deg)')
grid on
hold off

figure
hold on
plot (alpha*180/pi, '-b', 'LineWidth', 2)
xlabel ('Iteration')
ylabel ('Pitch angle (deg)')
grid on
hold off

figure
hold on
plot (y_L,theta*180/pi, '-b', 'LineWidth', 2)
xlabel ('Dimensionless span (y/L)')
ylabel ('Torsion angle (deg)')
grid on
hold off

figure
hold on
plot (y_L,wd,'-b', 'LineWidth', 2)
xlabel ('Dimensionless span (y/L)')
ylabel ('Bending (m)')
grid on
hold off