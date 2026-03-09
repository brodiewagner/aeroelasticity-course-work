% Simple VLM for general planforms, uses 1 chord-wise element
clear
clc

%% Warren-12 benchmark planform
c = 1.5 ;           % root chord
lambda = 2/3 ;      % taper
L = sqrt(2) ;       % semi span
LAM = 45 ;          % sweep angle of half chord [degrees]
S = 2*sqrt(2) ;     % use trapez rule for area if linear sweep
e = 0.25 ;          % elastic axis loc. rel. quarter chord

% inflow conditions
Vinf = 50 ;             % airspeed m/s
rho = 1.225 ;           % density kg/m^3
adeg = 6 ;              % root angle of attack degrees
alpha = adeg*pi/180 ;   % root angle of attack in radians

% numerical conditions
NP = 100 ;              % spanwise panels for vortex lattice, MUST BE EVEN!

% mesh the wing into NP panels and plot it out
DY = 2/(NP) ;                           % non dimensional panel span
Y_L = 0:DY:1 ;                          % non dimensional ordinate along one wing
y_L = -1:DY:1 ;                         % non dimensional ordinate along both wings (tip to tip)
C_Y = c*(1-lambda*Y_L) ;                % compute chord along span (linear taper only)
c_Y = [C_Y(length(C_Y):-1:2) C_Y] ;     % chord from tip to tip

% plot leading edge, trailing edge, 1/4 chord, 1/2 chord and 3/4 chord axes first half chord axis line from tip to tip:
x_hc = sign(y_L).*y_L*L*tand(LAM) ;             % half chord
x_le = -0.5*c_Y+x_hc ;                          % leading edge
x_te = 0.5*c_Y+x_hc ;                           % trailing edge
x_qc = -0.25*c_Y+x_hc ;                         % quarter chord
x_3qc = 0.25*c_Y+x_hc ;                         % 3/4 chord
z = zeros(size(x_hc)) ;                         % all points in a plane

%% compute points "A,B,C" for wing
A = [x_qc(1:NP); y_L(1:NP)*L; z(1:NP)] ;        % left vortex corner
B = [x_qc(2:NP+1); y_L(2:NP+1)*L; z(2:NP+1)] ;  % right vortex corner
Cy = 0.5*(A(2,:)+B(2,:)) ;                      % control points lie inbetween A and B for y coord
Cx = interp1(y_L,x_3qc,Cy/L) ;                  % interpolate x coordinate at 3/4 chord at these y points
C = [Cx; Cy; 0*Cy] ;

% vortex lattice method for wing consisting of NP panels with points and
% quarter chord A,B and collocation points at 3/4 chord C in incident flow
% V_inf at angle of attack alpha
 for j = 1:NP
    for k = 1:NP
 % call subroutine/function on to compute Velocity for unit vortex
 % strength at point C from vortex line from point A to point B
 % using Biot Savart Law
        VAB = V_AB(A(:,k),B(:,k),C(:,j)) ;
% call subroutine/function to compute Velocity for unit vortex
% strength at point C from vortex line from point at infinity to A
% using Biot Savart Law
        VAI = VA_INF(A(:,k),C(:,j)) ;
% call subroutine/function to compute Velocity for unit vortex
% strength at point C from vortex line from point B to infinity
% using Biot Savart Law
        VBI = VB_INF(B(:,k),C(:,j)) ;
% compute complete downwash velocity matrix for unit vortex strengths
        AIC(j,k,:) = VAB+VAI+VBI ;
    end
 end
 
AICx = squeeze(AIC(:,:,1)) ;    % get x component of downwash
AICy = squeeze(AIC(:,:,2)) ;    % get y component of downwash
AICz = squeeze(AIC(:,:,3)) ;    % get z component of downwash

% get unit normals for each panel
for k = 1:NP
    n(:,k) = cross(A(:,k)-C(:,k),A(:,k)-B(:,k)) ;   % unit normal from plane formed by A,B,C
    n(:,k) = n(:,k)/sqrt(dot(n(:,k),n(:,k))) ;
end
 
% form aerodynamic influence coefficient matrix projection of V on n
clear AIC
for i = 1:NP
    for j = 1:NP
        AIC(i,j) = n(1,j)*AICx(i,j) + n(2,j)*AICy(i,j) + n(3,j)*AICz(i,j) ;     % this is just a dot product
    end
end

% calculate incident flow velocity and apply V.n=0 boundary condition
V = -Vinf*(cos(alpha)*n(1,:)+sin(alpha)*n(3,:)) ;

% compute vortex strengths for given incident flow and wing geometry
gamma = AIC\V' ;

% these vortex strengths are at mid-points between A and B (vortex corners)
qinf = 0.5*rho*Vinf*Vinf ;      % dynamic head
Li = rho*Vinf*gamma.*DY*L ;     % calculate local lift
Ltot = sum(Li) ;                % total lift
Cltot = Ltot/(qinf*S) ;         % wing lift coefficient
Cla = Cltot/alpha ;             % Compute Lift-slope

% calculate local lift coefficients
for i = 1:NP
    localc(i) = (c_Y(i)+c_Y(i+1))/2. ;
    Cl(i) = Li(i)/(qinf*localc(i)*L*DY) ;
end

%% Plotting 

% wing geometry plotting
figure('Name', 'Warren-12 wing')
title('Warren-12 wing', 'FontSize', 14, 'FontWeight', 'bold')
plot3(x_le,y_L*L,z) ;                           % plot leading edge
hold on
plot3(x_hc,y_L*L,z,'--') ;                      % plot half chord
plot3(x_te,y_L*L,z) ;                           % plot trailing edge
plot3([x_le(NP+1) x_te(NP+1)],[L L],[0 0]) ;    % plot chord at tips
plot3([x_le(1) x_te(1)],[-L -L],[0 0]) ;        % plot chord at tips
plot3(x_qc,y_L*L,z,'r--')                       % plot out quarter chord, along which vortex positions lie
plot3(x_3qc,y_L*L,z,'g--')                      % plot out 3/4 chord, along which collocation points lie
grid on 
% plot points A,B,C as a check
plot3(A(1,:),A(2,:),A(3,:),'ro')
plot3(B(1,:),B(2,:),B(3,:),'r+')
plot3(C(1,:),C(2,:),C(3,:),'go')
hold off

% plot out CL/CL_total to show spanwise lift coefft distribution
figure('Name', 'Lift distribution for Warren-12 wing')
title('Span-wise Lift Coefficient for the Warren-12 wing', 'FontSize', 14, 'FontWeight', 'Bold')
hold on
grid on
plot (Cy((NP/2 + 1):NP),Cl((NP/2 + 1):NP)/Cltot, '-o', 'LineWidth', 2) ;
xlabel ('Span-wise position (-)', 'FontSize', 16) ;
ylabel ('Cl/Cl_{TOTAL}', 'FontSize', 16) ;

% plot out CL/CL_total to show spanwise lift coefft distribution
figure('Name', 'Loading for Warren-12 wing')
title('Span-wise Loading for the Warren-12 wing', 'FontSize', 14, 'FontWeight', 'Bold')
hold on
grid on
plot (Cy((NP/2 + 1):NP),Li((NP/2 + 1):NP), '-o', 'LineWidth', 2) ;
xlabel ('Span-wise position (-)', 'FontSize', 16) ;
ylabel ('Loading, $L_i$', 'FontSize', 16, 'Interpreter','latex') ;

