clc
clear

% Basic wing planform
c=1.5; %chord
lambda=0.4; % Taper ratio
L=6; % Semi span
e=0.25; % Elastic axis location relative to quarter chord
LAM=25; % sweep in degrees

% Meshing part
NP=100; % Number of panels
DY=2/(NP); % Dimensionless panel span
y_L=0:DY:1; % Panel position along one wing
c_Y=c*(1-lambda*y_L); % chord along the span due to taper
dy=DY*L; % Actual panel span
S=trapz(c_Y)*dy*2; % Wing area
AR=(L*2)^2/S; % Aspect ratio

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
Li=0.5*rho*Vinf^2*(c*(1-lambda*(y_L)))*a3*alpha(trimstep);
% Rigid wing CL along span
CL_y_R=Li./(0.5*rho*Vinf.^2*(c*(1-lambda*y_L)));

% step 2 Determine forces on each basis function via virtual work
Mi=e*c*Li;
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
Li=(0.5*rho*Vinf^2*(c*(1-lambda*(y_L)))*a3).*alphae';
Ltot=trapz(y_L,Li)*L*2; % calculate total wing lift
formatSpec = 'TrimStep %i, Ltot %f, mg %f, Lover %f, alpha %f, alphae(NP) %f\n';
fprintf (formatSpec,trimstep, Ltot, mg, Ltot-mg, alpha, alphae(NP/2))

% step 6 use simple trimming routine
for trimstep=2:200
% elastic angle of attack for one wing
alphae=alpha(trimstep-1)+theta*cosd(LAM)-wd*sind(LAM);
% lift for elastic wing
Li=(0.5*rho*Vinf^2*(c*(1-lambda*(y_L)))*a3).*alphae';
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
Mi=e*c*Li;
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
CL_y=Li./(0.5*rho*(Vinf^2)*(c*(1-lambda*(y_L))));
CL_sweep=(CL_y/CL)';

figure
hold on
plot (y_L,CL_y,'-b')
hold on
plot(y_L,CL_sweep,'-g')
hold on
plot (y_L,CL_y_R(1:51),'-r')
xlabel ('Dimensionless span (y/L)')
ylabel ('Lift coefficeint (CL)')
grid on
hold off

figure
hold on
plot (y_L,alphae*180/pi,'-b')
xlabel ('Dimensionless span (y/L)')
ylabel ('Aeroelastic angle (deg)')
grid on
hold off

figure
hold on
plot (alpha*180/pi,'-b')
xlabel ('Iteration')
ylabel ('Pitch angle (deg)')
grid on
hold off

figure
hold on
plot (y_L,theta*180/pi,'-b')
xlabel ('Dimensionless span (y/L)')
ylabel ('Torsion angle (deg)')
grid on
hold off

figure
hold on
plot (y_L,wd,'-b')
xlabel ('Dimensionless span (y/L)')
ylabel ('Bending (m)')
grid on
hold off