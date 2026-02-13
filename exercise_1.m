clc
clear all

%% example wing parameters
EA = 0.45 ;             % non dimensional elastic axis position chords
IA = 0.5 ;              % non dimensional inertial axis position chords
L = 6 ;                 % m wing semi span
c = 1 ;                 % m wing chord
b= c/2 ;                % m semi chord
xa = (IA-EA)*c/b ;      % wing chordwise cg semi-chords
mbar = 75 ;             % kg mass per unit span of wing
I = 25 ;                % kg m^2/m pitch moment of inertia of wing
EI = 2e6 ;              % Nm^2 flexural rigidity of wing
GJ = 5e5 ;              % Nm^2/rad torsional rigidity of wing

% store properties
mt = 100 ;              % kg mass of the store
Lt = 6 ;                % m spanwise position of store
It = 22.5 ;             % kg m^2 
x = 0.25 ;              % m position of the store ahead of the EA

% USING RAYLEIGH-RITZ FUNCTIONS OF w=a*(y/L)^2  AND theta=b*(y/L)
Mw = [mbar*L/5    0   ;
         0     I*L/3] ;
Mi = [mt    mt*x    ;
    mt*x It+mt*x^2] ;

% STIFFNESS MATRIX
K=[4*EI/L^3   0  ;
   0        GJ/L];
%echo off
% FREQUENCIES AND MODES OF WING+STORE FROM SYSTEM EIGENVALUES AND EIGENVECTORS
[V,LAMBDA] = eig(K,(Mw+Mi))
omega = sqrt(diag(LAMBDA));     % rad/s frequency
f = omega/(2*pi)                % hertz frequency
% echo off
% First two natural frequencies, Rayleigh-Ritz method, [Hz]:
% First mode=2.181 Hz
% Second mode=5.0969 Hz
% M-orthonormal modes, Rayleigh-Ritz method (Modes in columns, generalized coordinates in
% rows:
% V =
% 0.0696 -0.0352
% 0.0100 0.1084
% DEFLECTION AND TWIST OF ELASTIC AXIS:
 y_L=0:0.01:1;
 for i=1:2; % first two modes
     z_EA(i,:)=V(1,i)'*(y_L).^2;
 theta_EA(i,:)=V(2,i)'*(y_L);
     z_EA(i,:)=V(1,i)'*(y_L).^2;
 theta_EA(i,:)=V(2,i)'*(y_L);
 end;
figure 
plot(z_EA(1,:), 'LineWidth', 4);
hold on;
plot(z_EA(2,:), 'LineWidth', 4);
grid on 
grid minor

figure
plot(theta_EA(1,:)*180/pi, 'LineWidth', 4);
hold on;
plot(theta_EA(2,:)*180/pi, 'LineWidth', 4);
grid on 
grid minor